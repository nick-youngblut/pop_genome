#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $runID, $blast_params);
my $fork = 0;
my $num_threads = 1;
my $grouping = "-s";
GetOptions(
	   "fork=i" => \$fork,
	   "runID=s" => \$runID,
	   "blast=s" => \$blast_params,
	   "num_threads=i" => \$num_threads,
	   "grouping=s" => \$grouping,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a runID" unless $runID;


### MAIN
my $pm = new Parallel::ForkManager($fork);
while(<>){
	chomp;
	s/#.+//;
	next if /^\s*$/;
	
	# forking #
 	my $pid = $pm->start and next;
	
	# making blast query table from ITEP #
	my $query_file = genes_from_ITEP($_, $runID, $grouping);
	
	# blastp #
	call_blast($query_file, $blast_params, $num_threads);
	
	$pm->finish;
	}
$pm->wait_all_children;
	
### Subroutines
sub call_blast{
# calling blast on fasta #
	my ($query_file, $blast_params, $num_threads) = @_;
	
	(my $blast_out = $query_file) =~ s/\.[^\.]+$/.blast.txt/;
	$blast_params = "blastp -query $query_file -db nr -outfmt 5 -evalue 1e-5 -num_threads $num_threads | blastxml_parse.pl > $blast_out";
	print STDERR $blast_params, "\n" unless $verbose;
	
	`$blast_params`;
	}

sub genes_from_ITEP{
# getting genes from ITEP #
	my ($line, $runID, $grouping) = @_;
	
	# output file #
	my @line = split /\t/, $line;
	my $q = join("", join("\\t", @line), "\\n");
	$line[0] =~ s/ /_/g;
	my $outfile = "$line[0].faa";
	open OUT, ">$outfile" or die $!;

	# query ITEP #
	my $cmd = "printf \"$q\" | db_findClustersByOrganismList.py $runID $grouping | db_getGenesInClusters.py | db_getGeneInformation.py -g 3|";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, $cmd or die $!;
	while(<PIPE>){
		chomp;
		my @line = split /\t/;
		print OUT join("\n", ">$line[0]", $line[$#line]), "\n";
		}

	close PIPE;
	close OUT;
	
	return $outfile;
	}


__END__

=pod

=head1 NAME

blast_UniqueGenesFromITEP.pl -- blasting gene clusters only in 1 organim

=head1 SYNOPSIS

blast_UniqueGenesFromITEP.pl [flags] < organims 

=head2 required flags

=over

=item -runID

RunID for selecting genes (e.g. mazei_isolates_I_2.0_c_0.4_m_maxbit).

=back

=head2 optional flags

=over

=item -fork

Number of concurent blasts to run. [1]

=item -num_threads

Number of threads for each blast. [1]

=item blast

Blast parameters used. 
[blastp -query $query_file -db nr -outfmt 5 -evalue 1e-5]

-outfmt must be '5'

=item -grouping

Grouping parameter used for db_findClustersByOrganismList.py. [-s]

Default = genes only found in the organims provided.

See 'db_findClustersByOrganismList.py -h' for more info.

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc blast_GenesFromITEP.pl

=head1 DESCRIPTION

Blastp of all genes in unique gene clusters (clusters >= 1 in organism).

=head2 ITEP workflow for getting genes

db_findClustersByOrganismList.py --> db_getGenesInClusters.py --> db_getGeneInformation.py

=head2 Requires:

=over

=item * blastxml_parse.pl

=back

=head1 EXAMPLES

=head2 Usage (GTL isolates)

blast_genesFromITEP.pl <(grep "[1-3].[HF].[AMT]" organisms) -r mazei_isolates_I_2.0_c_0.4_m_maxbit -num 8

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

