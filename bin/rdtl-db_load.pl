#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
my ($db_file, $tree_meta_in, $node_in, $tree_in, $date);
my (@DTL_cost, $run_date, $include_leaves);
my $ranger_runID;
GetOptions(
	   "db=s" => \$db_file,
	   "metadata=s" => \$tree_meta_in, 
	   "runID=s" => \$ranger_runID,
	   "node=s" => \$node_in,
	   "tree=s" => \$tree_in,
	   "DTL_cost=i{3,3}" => \@DTL_cost,
	   "date=s" => \$run_date,
	   "leaf" => \$include_leaves,		# no
	   "replace" => \$replace,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
die " ERROR: provide a database file name!\n"
	unless $db_file;
die " ERROR: cannot find $db_file!\n"
	unless -e $db_file;


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
	or die " Can't connect to $db_file!\n";

load_metadata($dbh, $tree_meta_in) if $tree_meta_in;

if($node_in){			# loading node, tree, & runID, & costs
	# input check #
	my $runIDs_r = get_runIDs($dbh);
	check_runID($dbh, $runIDs_r, $ranger_runID);
	($tree_in, @DTL_cost, $run_date) = 
		check_input($tree_in, @DTL_cost, $run_date);	
	
	# loading ranger_run table #
	load_ranger_run($dbh, $ranger_runID, @DTL_cost, $run_date);
	
	# loading node table #
	load_node_table($dbh, $node_in, $ranger_runID);
	
	# loading tree table #
	load_tree_table($dbh, $tree_in, $ranger_runID);
	}

$dbh->commit;
$dbh->disconnect();
exit;


### Subroutines
sub load_tree_table{
# loading node table #
	my ($dbh, $tree_in, $ranger_runID) = @_;
	
	# sql #
	my $sql = "INSERT INTO Tree(Ranger_runID,TreeID,Min_rec_cost,D_cnt,T_cnt,L_cnt) values(?,?,?,?,?,?)";
	my $sth = $dbh->prepare($sql);
	
	# load table #
	open IN, $tree_in or die $!;
	my $entry_cnt = 0;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		$sth->execute($ranger_runID, @line);
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: loading line $. for tree file\n";
			}
		else{ $entry_cnt++; }
		}
	close IN;
	
	print STDERR "...Number of entries added to Tree table: $entry_cnt\n"
		unless $verbose;	
	}

sub load_node_table{
# loading node table #
	my ($dbh, $node_in, $ranger_runID) = @_;
	
	# sql #
	my $sql = "INSERT INTO Node(Ranger_runID,TreeID,Gene_nodeID,Species_nodeID,Category,T_recipient) values(?,?,?,?,?,?)";
	my $sth = $dbh->prepare($sql);
	
	# load table #
	open IN, $node_in or die $!;
	my $entry_cnt = 0;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		undef $line[4] unless $line[4];
		
		# skip leaves? #
		print Dumper $line[3]; exit;
		if($line[3] =~ /leaf/i){
			next unless $include_leaves;
			}
		
		# loading #
		$sth->execute(($ranger_runID, @line));
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: loading line $. for node file\n";
			}
		else{ $entry_cnt++; }
		}
	close IN;
	
	print STDERR "...Number of entries added to Node table: $entry_cnt\n"
		unless $verbose;	
	}

sub load_ranger_run{
# loading ranger_run table #
	my ($dbh, $ranger_runID, $D_cost, $T_cost, $L_cost, $run_date) = @_;
	
	my $sql = "INSERT INTO Ranger_run(Ranger_runID,D_cost,T_cost,L_cost,run_date) values(?,?,?,?,?)";
	my $sth = $dbh->prepare($sql);
	$sth->bind_param(1, $ranger_runID);
	$sth->bind_param(2, $D_cost);
	$sth->bind_param(3, $T_cost);
	$sth->bind_param(4, $L_cost);
	$sth->bind_param(5, $run_date);
	$sth->execute;
	if($DBI::err){
		print STDERR "ERROR: $DBI::errstr in: '", 
			join("\t", $ranger_runID, $D_cost, $T_cost, $L_cost, $run_date), 
			"'\n";
		}			
	}

sub load_metadata{
# loading metadata table #
	my ($dbh, $tree_meta_in) = @_;
	
	# sql #
	my $sql = "INSERT INTO Tree_meta(Cluster_runID,ClusterID,TreeID,Core_var,Date) VALUES(?,?,?,?,?)";
	my $sth = $dbh->prepare($sql);
	
	# reading in table; loading db #
	open IN, $tree_meta_in or die $!;
	my $entry_cnt = 0;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		$sth->execute(@line);
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '", join("\t", @line), "'\n";
			}
		else{ $entry_cnt++; }
		}
	close IN;

	print STDERR "...Number of entries added to Tree_meta table: $entry_cnt\n"
		unless $verbose;
	}

sub check_input{
# checking to make sure required input is provided #
	my ($tree_in, @DTL_cost, $run_date) = @_;
	
	# tree file #
	die " ERROR: provide a '*_tree.txt' file (from ranger-dtl_parse.pl)\n"
		unless $tree_in;
	
	# costs #
	die " ERROR: provide a duplication cost\n" unless $DTL_cost[0];
	die " ERROR: provide a transfer cost\n" unless $DTL_cost[1];
	die " ERROR: provide a loss cost\n" unless $DTL_cost[2];
	
	# date #
	unless($run_date){
		use POSIX qw(strftime);
		$run_date = strftime "%m/%d/%Y", localtime;
		}	
	
	return ($tree_in, @DTL_cost, $run_date);
	}

sub check_runID{
# checking for existence of runID #
	my ($dbh, $runIDs_r, $ranger_runID) = @_;
	die " ERROR: provide a Ranger-dtl run ID!\n"
		unless $ranger_runID;
		
	if(exists $runIDs_r->{$ranger_runID}){
		if($replace){
			
			my @tables = qw/ranger_runID node tree/;
			foreach my $table (@tables){
				my $sql = "DELETE FROM ? WHERE ranger_runID=?";
				my $sth = $dbh->prepare($sql);
				$sth->bind_param(1, $table);
				$sth->bind_param(2, $ranger_runID);
				$sth->execute;
				if($DBI::err){ print STDERR "ERROR: $DBI::errstr\n"; }
				}		

			print STDERR "...Ranger_runID: $ranger_runID deleted from rdtl-db!";
			}
		else{ 
			print STDERR " ERROR: runID -> '$ranger_runID' already found in DB!\n";
			print STDERR "### Existing RunIDs ###\n";
			print STDERR join("\n", keys @$runIDs_r), "\n";
			exit(1);
			}
		}
	}

sub get_runIDs{
# getting all unique runIDs from ranger_run table #
	my ($dbh) = @_;
	
	my $sql = "SELECT distinct(Ranger_runID) FROM Ranger_run";
	my $ret = $dbh->selectall_arrayref($sql);
	
	my %runIDs;
	map{$runIDs{$_} = 1 } @$ret;
	
	return \%runIDs;
	}


__END__

=pod

=head1 NAME

rdtl-db_load.pl -- load data into rdtl-db database

=head1 SYNOPSIS

rdtl-db_load.pl [flags]

=head2 Required flags

=over

=item -da 

rdtl-db database file.

=back

=head2 Optional flags

=over

=item -metadata

Tree metadata file (see discussion).

=item -node

'*_node.txt' from ranger-dtl_parse.pl

=item -tree

'*_tree.txt' from ranger-dtl_parse.pl

=item -runID

Ranger-DTL run ID. Any unique ID that you want.

=item -DTL_cost

3 arguments: duplication, transfer, loss cost used for Ranger-DTL run.

=item -date

Date of Ranger-DTL run. [today]

=item -replace

Replace Ranger-DTL run ID (and associated info)? [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc rdtl-db_load.pl

=head1 DESCRIPTION

Load rdtl-db with metadata on trees and/or data
from a Ranger-DTL run.

=head2 Columns needed in tree metadata file (no header in file)!

=over

=item *		ITEP Cluster run ID

=item * 	ITEP Cluster ID

=item * 	Tree ID (order in newick input to ranger-dtl)

=item * 	Core or variable gene tree (core|variable)

=item * 	Date (not required)

=back

=head1 EXAMPLES

=head2 Loading tree metadata:

rdtl-db_load.pl -db rdtl-db.sqlite -meta core_gene_trees.meta

=head2 Loading Ranger-DTL run data:

rdtl-db_load.pl -db rdtl-db.sqlite -run D1T1L1_core -DTL 1 1 1 -node ranger-dtl_parse_node.txt -tree ranger-dtl_parse_tree.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

