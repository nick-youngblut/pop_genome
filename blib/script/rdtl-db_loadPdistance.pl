#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Statistics::Descriptive;
use File::Temp qw/ tempfile tempdir /;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
my ($db_file, $align_list_in, $cluster_list_in);
my $move_id;
my $mothur_opts = "calc=onegap";
my $proc = 1;
GetOptions(
	   "database=s" => \$db_file,
	   "align=s" => \$align_list_in,
	   "cluster=s" => \$cluster_list_in,
	   "processors=i" => \$proc,
	   "mothur=s" => \$mothur_opts,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
check_file_IO($db_file, "database");
check_file_IO($align_list_in, "tree list");
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
my $aln_files_r = load_align_list($align_list_in);
my $clusters_r = load_cluster_list($cluster_list_in);
die " ERROR: number of alignment files != number of cluster names!\n"
	unless scalar @$aln_files_r == scalar @$clusters_r;

# calling mothur and getting pdistance stats #
my $tmpdir = File::Temp->newdir(); 		# temp directory
my $dirname = $tmpdir->dirname;
my $pdist_r = call_mothur_load_stats($dirname, $aln_files_r, $clusters_r);

# loading rdtl-db #
load_pdist_stats($dbh, $pdist_r);

# commit & disconnect #
$dbh->commit;
$dbh->disconnect();
exit;


### Subroutines
sub load_pdist_stats{
# loading bootstrap stats into rdtl-db #
	my ($dbh, $pdist_r) = @_;
	
	my @cols = qw/ClusterID Min Q1 Mean Median Q3 Max Stdev/;
	my $q = join("", 
			"INSERT INTO Pdistance(", join(",", @cols), 
			") VALUES(", join(",", ("?") x scalar @cols), ")"
			);
			
	my $sth = $dbh->prepare($q);
	
	my $entry_cnt = 0;
	foreach my $clusterID (keys %$pdist_r){
		$sth->execute($clusterID, @{$pdist_r->{$clusterID}});
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '$clusterID' \n";
			}
		else{ $entry_cnt++; }		
		}

	print STDERR " Number of entries added/updated in rdtl-db Pdistance table: $entry_cnt\n";		
	}

sub call_mothur_load_stats{
# calling mother (dist.seqs); getting stats from dist.seqs #
	my ($dirname, $aln_files_r, $clusters_r) = @_;
	
	my %pdist;
	for my $i (0..$#$aln_files_r){
		die " ERROR: cannot find cluster name for $$aln_files_r[$i]!\n"
			unless $$clusters_r[$i];
		
		# making symlink of file in temp directory #
		my $sl = make_tmp_symlink($$aln_files_r[$i], $dirname);
		
		# calling mothur #
		my $cmd = "mothur \"#dist.seqs(fasta=$sl, $mothur_opts, processors=$proc)\"";
		my $out = `$cmd`;
		print STDERR $out, "\n" if $verbose;
		
		# loading dist file; returning stats #		
		$pdist{$$clusters_r[$i]} = load_dist($sl);
		}

		#print Dumper %pdist; exit;
	return \%pdist;
	}
	
sub load_dist{
# loading distance, getting pdistances, calculating stats; returning stats #
	my $sl = shift;
	
	(my $dist = $sl) =~ s/\.[^.]+$|$/.dist/;
	die " ERROR: cannot find $dist!\n" unless -e $dist;
	
	my $stat = Statistics::Descriptive::Full->new();
	
	open IN, $dist or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split / /;
		die " ERROR: distance file does not have 3 columns!"
			unless scalar @l == 3;
			
		$stat->add_data($l[2]);
		}
	close IN;

	# getting stats #
	my @Q1 = $stat->percentile(25);
	my @Q3 = $stat->percentile(75);
	return [$stat->min(),
			$Q1[0],
			$stat->mean(),
			$stat->median(),
			$Q3[0],
			$stat->max(),
			$stat->standard_deviation() ];
	}
	
sub make_tmp_symlink{
	my ($file, $tmpdir) = @_;

	my @parts = File::Spec->splitpath($file);
	symlink($file, "$tmpdir/$parts[2]") or die $!;
	
	return "$tmpdir/$parts[2]";
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

sub load_align_list{
# loading a list of tree files #
	my ($align_list_in) = @_;
	
	open IN, $align_list_in or die $!;
	my %align_list;
	my @align_list;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		die " ERROR: $_ provided more than once!\n"
			if exists $align_list{$_};
		$align_list{$_} = 1;
		push @align_list, $_;
		}
	close IN;
	
	return \@align_list;
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

rdtl-db_loadPdistance.pl -- loading Pdistance stats trees of gene clusters

=head1 SYNOPSIS

rdtl-db_loadPdistance.pl [flags]

=head2 Required flags

=over

=item -database  <char>

rdtl-db database file.

=item -align  <char>

File with list of alignment files (1 per line; fasta format).

=item -cluster  <char>

File with list of clusters (1 per line; same order as '-align' list)

=back

=head2 Optional flags

=over

=item -mothur  <char>

Mothur dist.seqs paramters. ["calc=onegap"]

=item -processors  <int>

Number of processors used by mothur. [1]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc rdtl-db_loadPdistance.pl

=head1 DESCRIPTION

Load stats for Pdistance values of alignments
used for making trees, which were
used for the ranger-dtl run(s) in rdtl-db.

Provide a list of alignments and cluster IDs (same 
as in rdtl-db). The order of alignment files and 
clusterID must match!

=head1 EXAMPLES

=head2 Basic usage

rdtl-db_loadPdistance.pl -da rdtl-db.sqlite -t aln_file_list.txt -c clusterID_list.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

