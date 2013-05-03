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

my ($verbose, @sam_in, $database_file, $runID);
my $clusterID_col = 2;
GetOptions(
		"db=s" => \$database_file,		# database location
		"sam=s{,}" => \@sam_in,
		"column=i" => \$clusterID_col,
		"runID=s" => \$runID,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a run ID!\n" unless $runID;
map{ die " ERROR $_ not found!\n" unless -e $_ } @sam_in if @sam_in;

### MAIN
### Connect to DB
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!";

# loading #
my $contigs_r = load_contigs($dbh);
my $clusterID_r = load_cluster_ids($clusterID_col);
my $geneID_r = load_genes_from_clusters($clusterID_r, $runID);

foreach my $clusterID (keys %$clusterID_r){
	#my $gene_ids = get_geneIDs($clusterID, $runID);	
	
	#my $cluster_tbl = load_cluster_tbl($clusterID, $dbh);
	}

### Subroutines
sub load_cluster_tbl{
# loading cluster info needed for pulling out regions of interest #
	my ($clusterID, $dbh) = @_;
	
	#q1 = "SELECT geneid, genestart, geneend, strand, processed.contig_mod FROM processed WHERE geneid = $clusterID";
	

	}

sub load_genes_from_clusters{
# getting geneIDs for a cluster from ITEP #
	my ($clusterID_r, $runID) = @_;
	
	my @q;
	foreach my $clusterID (keys %$clusterID_r){
		push(@q, join("\\t", $runID, $clusterID));
		}
	my $q = join("\\n", @q);
	
	my $qq = "printf \"$q\" | db_getClusterGeneInformation.py |";
	print STDERR $q, "\n" unless $verbose;
	
	open PIPE, $qq or die $!;
	my %geneIDs;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		print Dumper @line;
		
		}
	
	close PIPE;
	exit;

	}

sub load_contigs{
# loading 'contig' table into memory #
	my ($dbh) = @_;
	
	my $q = "SELECT contig_mod, seq FROM contigs;";
	my $contig_seq_r = $dbh->selectall_arrayref($q); 
	
	my %contigs;
	foreach (@$contig_seq_r){
		$contigs{$$_[0]} = $$_[1];
		}
	
		#print Dumper %contigs; exit;
	return \%contigs;		# contig_mod => contig_seq
	}

sub load_cluster_ids{
# loading cluster_id file in 'ITEP' format #
	my ($clusterID_col) = @_;
	
	my %clusterID;
	while(<>){
		chomp;
		s/#.+//;
		next if /^\s$/;
		
		my @line = split /\t/;
		$clusterID{$line[$clusterID_col -1]} = 1;
		}
		
		print STDERR "Number of clusters provided: ", scalar keys %clusterID, "\n"
			unless $verbose;
			
		#print Dumper %clusterID; exit;
	return \%clusterID;
	}



__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

