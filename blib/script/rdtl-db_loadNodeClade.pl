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
my ($db_file, $taxa_clade_in);
my $delimit = "__";
GetOptions(
	   "database=s" => \$db_file,
	   "taxa=s" => \$taxa_clade_in,		
	   "delimiter=s" => \$delimit,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
die " ERROR: provide a database file name!\n"
	unless $db_file;
die " ERROR: cannot find $db_file!\n"
	unless -e $db_file;
die " ERROR: provide a taxa-clade list!\n"
	unless $taxa_clade_in;
die " ERROR: cannot find $taxa_clade_in!\n"
	unless -e $taxa_clade_in;


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
	or die " Can't connect to $db_file!\n";

# checking to make sure table exists #
my $table_list_r = list_tables($dbh);
die " ERROR: 'node_clade' table not found!\n"
	unless grep(/node_clade/i, @$table_list_r);

# loading taxa-clade list #
my $taxa_clade_r = load_taxa_clade_list($taxa_clade_in, $delimit);

# loading unique nodes #
my $sp_nodes_r = load_species_nodes($dbh);

# loading node_clade table #
load_node_clade_table($dbh, $taxa_clade_r, $sp_nodes_r);

$dbh->commit;
$dbh->disconnect();
exit;


### Subroutines
sub load_node_clade_table{
# loading node_clade table with unique nodes from nodes table #
	my ($dbh, $taxa_clade_r, $sp_nodes_r) = @_;
	
		#print Dumper %$sp_nodes_r; exit;
	
	# sql #
	my $sql = "INSERT INTO node_clade(Species_NodeID, CladeID) values(?,?)";
	my $sth = $dbh->prepare($sql);
	
	my $entry_cnt = 0;
	foreach my $sp_node (keys %$sp_nodes_r){
		my @taxa = split /\|/, $sp_node;
		die " ERROR: number of nodes != 1 or 2 for $sp_node!\n"
			unless scalar @taxa == 1 || scalar @taxa == 2;
		
		my %clades;
		my $next = 0;
		foreach (@taxa){
			unless (exists $taxa_clade_r->{$_}){
				print STDERR " WARNING: cannot find $_ in taxa_clade list! Skipping!\n";
				$next++;
				next;
				}			
			$clades{$taxa_clade_r->{$_}} = 1;
			}
		next if $next;
		
			#print Dumper $sp_node;
			#$sth->bind_param(1, $sp_node);
		
		if(scalar keys %clades == 1){		# node is monophyletic for specified clades
			my @clades = keys %clades;
			#$sth->bind_param(2, $clades[0]);
			$sth->execute($sp_node, $clades[0]);
			}
		else{				# node is polyphyletic 
			$sth->execute($sp_node, my $tmp);
			}
			

		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '$sp_node' \n";
			}
		else{ $entry_cnt++; }
		}
	
	print STDERR " Number of entries added/updated in rdtl-db node_clade table: $entry_cnt\n";
	}

sub load_species_nodes{
# loading all unique species nodes from the Nodes table in rdtl-db #
	my ($dbh) = @_;

	my $q = "SELECT Species_NodeID, T_recipient FROM Node GROUP BY Species_NodeID, T_recipient";
	my $ret = $dbh->selectall_arrayref($q);
	die " ERROR: no matching entries!\n" unless @$ret;
	
	my %sp_nodes;
	foreach my $entry (@$ret){
		foreach (@$entry){
			next unless $_;			# if not T_recipent
			$sp_nodes{$_} = 1;
			}
		}
		
		#print Dumper %sp_nodes; exit;
	return \%sp_nodes;
	}

sub load_taxa_clade_list{
# loading taxa clade list #
# taxa'delimit'clade #
	my ($taxa_clade_in, $delimit) = @_;
	
	open IN, $taxa_clade_in or die $!;
	my %taxa_clade;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /$delimit/;
		die " ERROR: '$_' is not split into 2 columns by '$delimit'!\n"
			unless scalar @l == 2;

		$taxa_clade{$l[0]} = $l[1]; 
		}
	close IN;
	
		#print Dumper %taxa_clade; exit;
	return \%taxa_clade;
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

sub list_tables{
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}

__END__

=pod

=head1 NAME

rdtl-db_loadNodeClade.pl -- loading clade IDs for nodes in rdtl-db

=head1 SYNOPSIS

rdtl-db_loadNodeClade.pl [flags]

=head2 Required flags

=over

=item -database  <char>

rdtl-db database file.

=item -taxa  <char>

taxa-clade list file (tab-delimited; 2 column: taxon, clade).

=back

=head2 Optional flags

=over

=item -delimiter  <char>

Delimiter for parsing taxa-clade list into 2 columns. ['__']

=item -h	This help message

=back

=head2 For more information:

perldoc rdtl-db_loadNodeClade.pl

=head1 DESCRIPTION

Load cladeIDs for nodes (both internal & external)
into rdtl-db. Clade IDs are determine from a list
of extant taxa (tree leaves) and which clade they 
fall into. rename_addCladeName.pl can be used to make
the list. Internal nodes are classified by clade
based on their extant child nodes. Nodes polyphyletic
to the clades specified will have NULL cladeIDs.

The node-clade list should be based off of the 
species tree used for the ranger-dtl runs.

=head1 EXAMPLES

=head2 Method for making node-clade list table

rename_addCladeName.pl <(nw_labels -I species.nwk) -t species.nwk -r 1 10 11 -n clade1 clade2 > taxa_clade.txt

=head2 Loading database

rdtl-db_makeNodeClade.pl -da rdtl-db.sqlite -t taxa_clade.txt 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

