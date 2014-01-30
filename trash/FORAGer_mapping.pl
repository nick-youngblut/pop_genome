#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $qdir, $sdir, $regex, $list_in);
my $bowtie_params = "-k 10";
my $fork = 0;
GetOptions(
	   "query=s" => \$qdir,			# query directory
	   "subject=s" => \$sdir,		# subject directory
	   "list=s" => \$list_in, 		# list of queries
	   "forks=i" => \$fork, 		# number of forked mappings
	   "params=s" => \$bowtie_params,	# params passed to bowtie2
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a query (read files) directory!\n" unless $qdir;
die " ERROR: provide a subject (reference genomes) directory!\n" unless $sdir;
map{ $_ = File::Spec->rel2abs($_) } ($qdir, $sdir);
map{ die " ERROR: $_ not found!\n" unless -d $_ } ($qdir, $sdir);
die " ERROR: provide a list file!\n" unless $list_in;
die " ERROR: $list_in not found!\n" unless -e $list_in;

### MAIN
# loading directory files #
my $qfiles_r = load_dir($qdir);
my $sfiles_r = load_dir($sdir);

# loading list #
my $qlist_r = load_list($list_in);


# building indices for bowtie2 #
my $pmb = new Parallel::ForkManager($fork);
foreach my $subject (@$sfiles_r){
	my $pid = $pmb->start and next;
	call_bowtie2_build($subject, $sdir);
	$pmb->finish;
	}
$pmb->wait_all_children;

# making output dir #
my $outdir = make_outdir();

# read mapping #
my $pm = new Parallel::ForkManager($fork);
foreach my $org (keys %$qlist_r){			# each set of read files
	foreach my $subject (@$sfiles_r){		# each subject genome
		
		#forking #
		my $pid = $pm->start and next;
	
		# bowtie2 mapping #
		call_bowtie2($org, $subject, $qlist_r, $bowtie_params, 
					$sdir, $qdir, $outdir);
	
		# end fork #
		$pm->finish;
		}
	}
$pm->wait_all_children;



### Subroutines
sub make_outdir{
	my $outdir = File::Spec->rel2abs("./sam/");
	rmtree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	
	print STDERR "Writing files to: $outdir\n";
	return $outdir;
	}

sub call_bowtie2{
	my ($org, $subject, $qlist_r, $bowtie_params, 
					$sdir, $qdir, $outdir) = @_;
	
	my $cmd;
	if(exists $qlist_r->{$org}{"F"} && exists $qlist_r->{$org}{"R"}){	# 2 read files
		$cmd = "bowtie2 -x $sdir/$subject -S $outdir/$org-$subject.sam -1 $qdir/$qlist_r->{$org}{'F'} -2 $qdir/$qlist_r->{$org}{'R'} $bowtie_params"; 
		}
	elsif(exists $qlist_r->{$org}{"FR"}){
		$cmd = "bowtie2 -x $sdir/$subject -S $outdir/$org-$subject.sam -U $qdir/$qlist_r->{$org}{'FR'} $bowtie_params"; 
		}
	print STDERR $cmd, "\n";
	`$cmd`;
	}

sub call_bowtie2_build{
	my ($subject, $sdir) = @_;
	
	my $cmd = "bowtie2-build $sdir/$subject $sdir/$subject";
	print STDERR $cmd, "\n";
	`$cmd`;
	}	

sub load_list{
	my ($list_in) = @_;
	
	open IN, $list_in or die $!;
	my %qlist;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: $line[0] needs 1 or 2 reads file names in the list file!\n"
			unless scalar @line >=2;
		if(scalar @line == 2){ $qlist{$line[0]}{"FR"} = $line[1]; }
		elsif(scalar @line == 3){
			$qlist{$line[0]}{"F"} = $line[1];
			$qlist{$line[0]}{"R"} = $line[2];
			}
		else{ die " ERROR: $!\n"; }
		}
	close IN;
	
	return \%qlist;
	}

sub load_dir{
	my ($dir) = @_;
	opendir DIN, $dir or die $!;
	my @files = grep(/\.(fq|fastq|fasta|fna|fa)$/, readdir DIN);
	closedir DIN;
	
	return \@files;
	}



__END__

=pod

=head1 NAME

FORAGer_mapping.pl -- Mapping reads from query genomes to all reference genomes

=head1 SYNOPSIS

FORAGer_mapping.pl [flags]

=head2 Required flags

=over

=item -query

Directory of query read files (.fq|.fastq|.fasta|.fna|.fa).
Use '-f' for bowtie2 params if files are fasta format.

=item -subject

Directory of reference genomes (fasta format; extension = '.fasta|.fna|.fa').

=item -list

List file; tab-delimited; no header

3 columns: query_organism_name, read1_file_name, read2_file_name

=back

=head2 Optional flags

=over

=item -params

Parameters passed to bowtie2 (besides '-x, -S, -1, -2'). [-k 10]

=item -forks

Number of parallel calls of bowtie2. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_mapping.pl

=head1 DESCRIPTION

Just a wrapper for bowtie2 to facilitate mapping of the query reads
to all of the subject (reference) genomes.

Output SAM file naming: 'query'-'subject'.sam

By default, the top 10 hits for each read are kept (bowtie2 param: '-k 10'). 

=head1 EXAMPLES

=head2 Basic Usage

FORAGer_mapping.pl -q ./query/ -s ./subject/ -list query_list.txt

=head2 Forking & multiple bowtie2 threads

FORAGer_mapping.pl -q ./query/ -s ./subject/ -list query_list.txt -f 5 -p "-k 10 -p 4"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

