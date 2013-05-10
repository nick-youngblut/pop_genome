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

my ($verbose, $nuc_bool, $blast_params, $lower_tri, $diag);
my $fork = 0;
my $cutoff = 100;
GetOptions(
	   "extra=s" => \$blast_params,		# blast parameters
	   "nucleotide" => \$nuc_bool, 		# nucleotide? [false]
	   "cutoff=f" => \$cutoff, 			# percent similarity to be considered the same
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
# loading fasta #
chomp(@ARGV = <>) unless @ARGV;

# running blast #
my $pm = new Parallel::ForkManager($fork);
foreach my $fasta_in (@ARGV){
	$pm->start and next;
	
	die " ERROR: $fasta_in not found!\n" unless -e $fasta_in;

	my $blast_res_r = parse_blast($fasta_in, $nuc_bool, $blast_params);	
	write_blast_sum_tbl($blast_res_r, $cutoff, $fasta_in);
	
	$pm->finish;
	}
$pm->wait_all_children;

### Subroutines
sub write_blast_sum_tbl{
	my ($blast_res_r, $cutoff, $fasta_in) = @_;
	
	my @keys = keys %$blast_res_r;
	my ($N_invar, $N_total) = (0,0);
	foreach my $i (0..$#keys){
		foreach my $ii (0..$#keys){
			next if $ii >= $i;			# no diag
			unless (exists $blast_res_r->{$keys[$i]}{$keys[$ii]}){
				print STDERR " WARNING: blast comparison $keys[$i] <-> $keys[$ii] not found!\n";
				next;
				}
			
			$N_invar++ if $blast_res_r->{$keys[$i]}{$keys[$ii]} >= $cutoff;
			$N_total++;
			}
		}
	print join("\t", $fasta_in, $N_invar, $N_total, sprintf("%.2f", $N_invar/$N_total * 100)), "\n";
	}

sub parse_blast{
	my ($fasta_in, $nuc_bool, $blast_params) = @_;
	my $blast;
	$nuc_bool ? $blast = "blastn" : $blast = "blastp";
	open PIPE, "$blast -subject $fasta_in -query $fasta_in -outfmt 6 |";
	
	my %blast_res;
	while(<PIPE>){
		chomp;
		my @line = split /\t/;
		$blast_res{$line[0]}{$line[1]} = $line[2];
		}
	close PIPE;
	
		#print Dumper %blast_res; exit;
	return \%blast_res;
	}


__END__

=pod

=head1 NAME

geneClusterInvariant.pl -- determine the percentage of genes in a cluster that are the same

=head1 SYNOPSIS

geneClusterInvariant.pl [options] < *fasta > invar.txt

=head2 options

=over

=item -cutoff <float>

The sequence identity cutoff for calling two sequences the 'same'. [100]

=item -nucleotide <bool>

blastn or blastp. [blastp]

=item -extra <char>

Extra blast parameters. []

=item -fork <integer>

Number of parallel blast runs. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc geneClusterInvariant.pl

=head1 DESCRIPTION

Use blast to determine how many genes in a fasta (nuc | AA) are 
the 'same'. By default 'same' = 100% sequence identity for 
what blast. 

Blast parameters used: "blast(p|n) -subject -query -outfmt 6"

Output: file N_sequences N_'same' Percent_'same'

=head1 EXAMPLES

=head2 Usage (4 parallel blastp)

find . -name "*fasta" | geneClusterInvariant.pl -f 4

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

