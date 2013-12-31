#!/usr/bin/env perl

=pod

=head1 NAME

PopGen-db_loadClusterInfo.pl -- loading cluster info from an ITEP geneInfo table

=head1 SYNOPSIS

PopGen-db_loadClusterInfo.pl [flags] < geneInfo.txt

=head2 Required flags

=over

=item -database  <char>

PopGen-db database file.

=item -runID  <char>

ITEP cluster run ID. 

=item -core  <int>

ClusterID considered core if found in '-core' taxa.

=back

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc PopGen-db_loadClusterInfo.pl

=head1 DESCRIPTION

Load basic cluster info into PopGen-db.
The table produced by db_getClusterGeneInformation.py 
is used for loading the database.

=head1 EXAMPLES

=head2 Basic usage

PopGen-db_loadClusterInfo.pl -da PopGen-db.sqlite -r mazei_I_2.0_c_0.4_m_maxbit -c 56 < gene_info.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


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

my ($verbose, $db_file, $core, $runID);
GetOptions(
	   "database=s" => \$db_file,
	   "core=i" => \$core,
	   "runID=s" => \$runID,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
check_file_IO($db_file, "database");
die " ERROR: provide a value for '-core'\n" unless defined $core;
die " ERROR: provide a value for '-runID'\n" unless defined $runID;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
	or die " Can't connect to $db_file!\n";

# loading gene info #
my $info_r = load_geneInfo();

# loading PopGen-db #
load_clust_info($dbh, $info_r, $runID, $core);

# commit & disconnect #
$dbh->commit;
$dbh->disconnect();
exit;


### Subroutines
sub load_clust_info{
# loading bootstrap stats into rdtl-db #
	my ($dbh, $info_r, $runID, $core) = @_;
	
	my @cols = qw/clusterID runID core_var/;
	my $q = join("", 
			"INSERT INTO Cluster_meta(", join(",", @cols), 
			") VALUES(", join(",", ("?") x scalar @cols), ")"
			);
			
	my $sth = $dbh->prepare($q);
	
	my %entry_cnt;
	foreach my $clusterID (keys %$info_r){
		
		# determining core #
		my $core_var = "variable";
		my $Nuniq = 0;
		map{$Nuniq++ if $info_r->{$clusterID}{$_} == 1} keys %{$info_r->{$clusterID}};
		$core_var = "core" if $Nuniq == $core;
		
		$sth->execute($clusterID, $runID, $core_var);
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '$clusterID' \n";
			}
		else{
			 $entry_cnt{$core_var}++; 
			 $entry_cnt{'total'}++; 
			 }		
		}

	map{ $entry_cnt{$_} = 0 unless exists $entry_cnt{$_} } qw/core variable total/;
	print STDERR " Number of entries added/updated in PopGen-db Cluster_meta table: ",
					$entry_cnt{'total'}, "\n";		
	print STDERR " \tNumber of 'core' entries: ",
					$entry_cnt{'core'}, "\n";		
	print STDERR " \tNumber of 'variable' entries: ",
					$entry_cnt{'variable'}, "\n";		

	}

sub load_geneInfo{
	my %info;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		die " ERROR: no clusterID column (column 14) in geneInfo table!\n"
			unless $l[13] && $l[13] =~ /^\d+$/;
		$info{$l[13]}{$l[2]}++;
		}	
		#print Dumper %info; exit;
	return \%info;
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



