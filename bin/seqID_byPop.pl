#!/usr/bin/env perl

=pod

=head1 NAME

seqID_byPop.pl -- Calculating average sequence identity within & between populations

=head1 SYNOPSIS

seqID_byPop.pl [flags] > seqID_stats.txt

=head2 Required flags

=over

=item -table  <char>

Table of tree files & matching gene cluster IDs (2-column; tab-delimited; 'file\tclusterID').

=item -population  <char>

Table of taxa in tree files & matching population (2-column; tab-delimited; 'taxon\tpopulation')

=back

=head2 Optional flags

=over

=item -delimiter  <char>

Delimiter separating taxon_name from annotation in sequence files. [" "]

=item -processors  <int>

Number of processors used by Mothur. [1]

=item -mothur  <char>

Mothur cmd used to produce a pdistance matrix. ["dist.seqs(fasta=?, calc=onegap, countends=T)"]

=item -h	This help message

=back

=head2 For more information:

perldoc seqID_byPop.pl

=head1 DESCRIPTION

Get average percent sequence identity of all comparisons of sequences.
Average percent identity comparisons are parsed by population
(both within and between) if a population file is provided. 

Not all taxa in the population file need to be found in each alignment.
The '-delimiter' option can be used is annotations or other info is 
appended on the taxon names in the sequence files.

=head1 EXAMPLES

=head2 Basic usage

seqID_byPop.pl -p pop.txt -t file-cluster_list.txt > seqID_stats.txt

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
use Bio::AlignIO;
use Statistics::Descriptive;
use File::Path qw/remove_tree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $pop_in, $file_in);
my $procs = 1;
my $delim = " ";
my $mothur_cmd = "dist.seqs(fasta=?, calc=onegap, countends=T)";
GetOptions(
	   "table=s" => \$file_in,
	   "population=s" => \$pop_in, 	# population table file
	   "delimiter=s" => \$delim,
	   "processors=i" => \$procs,
	   "mothur=s" => \$mothur_cmd,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$delim = qr/$delim/;
die "ERROR: provide a population file (-p)!\n" unless defined $pop_in;
die "ERROR: provide a file table (-t)!\n" unless defined $file_in;
die "ERROR: cannot find $file_in" unless -e $file_in;


### MAIN
# getting working directory #
my $cwd = File::Spec->rel2abs(File::Spec->curdir());

# loading file table #
my $files_r = load_file_table($file_in);

# loading population file (if provided) #
my $pops_r = load_pop($pop_in);

# processing each input file #
print join("\t", qw/file cluster population min q1 mean median q3 max N/), "\n";
foreach my $infile (keys %$files_r){
	$infile = File::Spec->rel2abs($infile);
		
	# making temp directory #
	use File::Temp qw/tempdir/;
	my $tmpdir = File::Temp->newdir(); 		# temp directory
	my $dirname = $tmpdir->dirname;
	chdir($dirname) or die $!;
	
	# symlinking #
	my @parts = File::Spec->splitpath($infile);
	my $link_name = "$dirname/$parts[2]";
	symlink($infile, $link_name) or die $!;
	
	# calling mother #
	my $dist_file = call_mother($link_name, $mothur_cmd, $procs);
	
	# loading distance file & summing distances per population comparison #
	my $stats_r = sum_distances($dist_file, $pops_r, $infile);
	
	# writing out table #
	write_stats_table($stats_r, $files_r);
	
	# moving back to original cwd #
	chdir($cwd) or die $!;	
	}


### Subroutines
sub write_stats_table{
	my ($stats_r, $files_r) = @_;
	
	my @cat = qw/min q1 mean median q3 max N/;
	#print join("\t", qw/file cluster population/, @cat), "\n";
	foreach my $file (keys %$stats_r){
		foreach my $pop (keys %{$stats_r->{$file}}){
			# NA unless value for each stat category #
			map{ $stats_r->{$file}{$pop}{$_} = "NA" unless
					exists $stats_r->{$file}{$pop}{$_} } @cat;
					
			print join("\t", $file, $files_r->{$file}, $pop, @{$stats_r->{$file}{$pop}}{@cat}), "\n";
			}
		}
	}

sub sum_distances{
	my ($dist_file, $pops_r, $infile) = @_;
	
	# getting unique pops #
	my %upop;
	map{ $upop{$_} = 1 } values %$pops_r;
	my @upop = sort keys %upop;
	my @pop_comb;
	for my $i (0..$#upop){
		for my $ii (0..$#upop){
			next if $ii < $i;
			push @pop_comb, join("__", $upop[$i], $upop[$ii]);
			}
		}
	push @pop_comb, "total";

	# loading distances from Mothur file and summing (to get mean) #
	my %pop_seqID;
	open IN, $dist_file or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split / /;
		
		for my $i (0..1){
			my @tmp = split /$delim/, $line[$i], 2;
			$line[$i] = $tmp[0];
			}
		
		# pdist to seqID #
		$line[2] = 100 - $line[2] * 100;
		
		# checking for pop of each taxon #
		die " ERROR: $dist_file -> $line[0] not found in population file!\n"
			unless exists $pops_r->{$line[0]};
		die " ERROR: $dist_file -> $line[1] not found in population file!\n"
			unless exists $pops_r->{$line[1]};

		# loading seqID into correct population comparison #
		my $pops = join("__", sort{$a cmp $b} ($pops_r->{$line[0]}, $pops_r->{$line[1]}) );


		# initializing stat objects #
		$pop_seqID{$pops} = Statistics::Descriptive::Full->new()
			unless exists $pop_seqID{$pops};				
		$pop_seqID{'total'} = Statistics::Descriptive::Full->new()
			unless exists $pop_seqID{'total'};
		
		# adding data
		$pop_seqID{$pops}->add_data($line[2]);
		$pop_seqID{'total'}->add_data($line[2]);
		
		}
	close IN;
	
	# seqID stats #
	my %stats;
	foreach my $pop (keys %pop_seqID){
		my ($min, $q1, $mean, $median, $q3, $max, $N) = ('NA') x 7;
		$min = $pop_seqID{$pop}->min() if defined $pop_seqID{$pop}->min();
		$q1 = ($pop_seqID{$pop}->percentile(25))[0] if defined $pop_seqID{$pop}->percentile(25);
		$mean = $pop_seqID{$pop}->mean() if defined $pop_seqID{$pop}->mean();
		$median = $pop_seqID{$pop}->median() if defined $pop_seqID{$pop}->median();
		$q3 = ($pop_seqID{$pop}->percentile(75))[0] if defined $pop_seqID{$pop}->percentile(75);
		$max = $pop_seqID{$pop}->max() if defined $pop_seqID{$pop}->max();
		$N = $pop_seqID{$pop}->count() if defined $pop_seqID{$pop}->count();
		
		$stats{$infile}{$pop} = {
			min => $min,
			q1 => $q1,
			mean => $mean,
			median => $median,
			q3 => $q3,
			max => $max,
			N => $N
			};
		}

		#print Dumper %stats; exit;
	return \%stats;
	}

sub load_file_table{
# loading file table
# 	2 columns: file\tcluster
	my ($file_in) = @_;
	
	open IN, $file_in or die $!;
	my %files;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die "ERROR: '$_' not in 2-column format!\n"
			unless scalar @l == 2;

		$l[0] = File::Spec->rel2abs($l[0]);
		$files{$l[0]} = $l[1];
		}
	close IN;
	
		#print %files;
	return \%files;
	}

sub call_mother{
# calling mothur to get pairwise distances among all sequences #
	my ($infile, $mothur_cmd, $procs) = @_;
	$mothur_cmd =~ s/\?/$infile/g;
	$mothur_cmd =~ s/\)/, processors=$procs)/;
		#print Dumper $mothur_cmd; exit;
	`mothur "#$mothur_cmd"`;
	(my $dist_file = $infile) =~ s/\.[^.]+$|$/.dist/;
	return $dist_file;
	}

sub load_pop{
# loading population table #
	my ($pop_in) = @_;
	open IN, $pop_in or die $!;
	my %pops;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: population table is not formatted correctly\n"
			unless scalar @line == 2;
		$pops{$line[0]} = $line[1];
		}
	close IN;
		#print Dumper %pops; exit;
	return \%pops;
	}
	

