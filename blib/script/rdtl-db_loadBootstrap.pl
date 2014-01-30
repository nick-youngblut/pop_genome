#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Bio::TreeIO;
use Statistics::Descriptive;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
my ($db_file, $tree_list_in, $cluster_list_in);
my $move_id;
GetOptions(
	   "database=s" => \$db_file,
	   "tree=s" => \$tree_list_in,
	   "cluster=s" => \$cluster_list_in,
	   "id" => \$move_id,			# move node ID to bootstrap? [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
check_file_IO($db_file, "database");
check_file_IO($tree_list_in, "tree list");
check_file_IO($cluster_list_in, "cluster list");

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
	or die " Can't connect to $db_file!\n";

# checking to make sure table exists #
my $table_list_r = list_tables($dbh);
die " ERROR: 'bootstrap' table not found!\n"
	unless grep(/bootstrap/i, @$table_list_r);
	
# loading lists #
my $tree_files_r = load_tree_list($tree_list_in);
my $clusters_r = load_cluster_list($cluster_list_in);
die " ERROR: number of tree files != number of cluster names!\n"
	unless scalar @$tree_files_r == scalar @$clusters_r;

# loading trees & calculating bootstrap stats #
my $boot_r = get_bootstrap_stats($tree_files_r, $clusters_r);

# loading stats into rdtl-db #
load_boot_stats($dbh, $boot_r);


# commit & disconnect #
$dbh->commit;
$dbh->disconnect();
exit;


### Subroutines
sub load_boot_stats{
# loading bootstrap stats into rdtl-db #
	my ($dbh, $boot_r) = @_;
	
	my @cols = qw/ClusterID Min Q1 Mean Median Q3 Max Stdev/;
	my $q = join("", 
			"INSERT INTO Bootstrap(", join(",", @cols), 
			") VALUES(", join(",", ("?") x scalar @cols), ")"
			);
			
	my $sth = $dbh->prepare($q);
	
	my $entry_cnt = 0;
	foreach my $clusterID (keys %$boot_r){
		$sth->execute($clusterID, @{$boot_r->{$clusterID}});
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '$clusterID' \n";
			}
		else{ $entry_cnt++; }		
		}

	print STDERR " Number of entries added/updated in rdtl-db bootstrap table: $entry_cnt\n";		
	}

sub get_bootstrap_stats{
# loading bootstraps for each tree and calculating stats #
	my ($tree_files_r, $clusters_r) = @_;
	
	my %boot;
	for my $i (0..$#$tree_files_r){
		
		# status #
		print STDERR "...processing: $$tree_files_r[$i]\n"
			if $verbose;
		
		# sanity check #
		die " ERROR: cannot find cluster for $$tree_files_r[$i]!\n"
			unless $$clusters_r[$i];
		
		# loading tree #
		my $treeo = tree_io($$tree_files_r[$i]);
		$treeo->move_id_to_bootstrap unless $move_id;
		
		# loading bootstraps into stats object #
		my $stat = Statistics::Descriptive::Full->new();
		for my $node ($treeo->get_nodes){
			next if $node->is_Leaf;
			$node->bootstrap(0) unless $node->bootstrap;
			
				#print Dumper join("__", $node->id, $node->bootstrap);
			die " ERROR: no bootstrap support found to a node in $$tree_files_r[$i]!\n"
				unless $node->bootstrap =~ /^\d+$/;
			$stat->add_data($node->bootstrap);
			}
			
		# getting stats #
		my @Q1 = $stat->percentile(25);
		my @Q3 = $stat->percentile(75);
		$boot{$$clusters_r[$i]} = [
					$stat->min(),
					$Q1[0],
					$stat->mean(),
					$stat->median(),
					$Q3[0],
					$stat->max(),
					$stat->standard_deviation() ];
		
		}
		#print Dumper %boot; exit;
	return \%boot;
	}

sub tree_io{
	# loading tree object: just 1st tree #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => "newick");
	my $treeio = $input->next_tree;	
	return $treeio;
	}

sub load_cluster_list{
# loading list of gene clsuters corresponding to the tree cluster list #
	my ($cluster_list_in) = @_;
	
	open IN, $cluster_list_in or die $!;
	my %clusters;
	my @clusters;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		die " ERROR: $_ provided more than once!\n"
			if exists $clusters{$_};
		$clusters{$_} = 1;
		push @clusters, $_;
		}
	close IN;
	
	return \@clusters;
	}

sub load_tree_list{
# loading a list of tree files #
	my ($tree_list_in) = @_;
	
	open IN, $tree_list_in or die $!;
	my %tree_list;
	my @tree_list;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		die " ERROR: $_ provided more than once!\n"
			if exists $tree_list{$_};
		$tree_list{$_} = 1;
		push @tree_list, $_;
		}
	close IN;
	
	return \@tree_list;
	}

sub check_file_IO{
	my ($infile, $var) = @_;
	die " ERROR: provide a $var file!\n"
		unless $infile;
	die " ERROR: cannot find $infile!\n"
		unless -e $infile;
	}

sub list_tables{
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}

__END__

=pod

=head1 NAME

rdtl-db_loadBootstrap.pl -- loading bootstrap stats trees of gene clusters

=head1 SYNOPSIS

rdtl-db_loadBootstrap.pl [flags]

=head2 Required flags

=over

=item -database  <char>

rdtl-db database file.

=item -tree  <char>

File with list of tree files (1 per line; newick files).

=item -cluster  <char>

File with list of clusters (1 per line; same order as '-tree' list)

=back

=head2 Optional flags

=over

=item -id  <bool>

Are bootstrap values encoded as internal node IDs? [TRUE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc rdtl-db_loadBootstrap.pl

=head1 DESCRIPTION

Load stats for boostrap values of trees 
used for the ranger-dtl run(s) in rdtl-db.

Provide a list of trees and cluster IDs (same 
as in rdtl-db). The order of tree files and 
clusterID must match!

=head1 EXAMPLES

=head2 Basic usage

rdtl-db_loadBootstrap.pl -da rdtl-db.sqlite -t tree_file_list.txt -c clusterID_list.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

