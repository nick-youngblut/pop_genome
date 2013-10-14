#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempfile tempdir /;
use Forks::Super;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $all_tests;
my $phi_window = 100;
my $phi_perm = 1000;
my $fork = 1;
GetOptions(
	   "o" => \$all_tests,			# TRUE
	   "w" => \$phi_window,
	   "p" => \$phi_perm,
	   "fork=i" => \$fork, 		
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide >=1 fasta file name!\n"
	unless $ARGV[0];
map{" ERROR: cannot find $_!\n" unless -e $_ } @ARGV;

### MAIN
# output header #
print join("\t", qw/file start end test pvalue/), "\n";

my %res_all;
foreach my $infile (@ARGV){
	my $job = fork {
			max_proc => $fork,
			share => [ \%res_all ],
			sub => sub {
				
				# calling phi #
				my $cmd = "Phi -f $infile -w $phi_window";
				$cmd .= " -p $phi_perm" if $phi_perm;
				$cmd .= " -o" if $all_tests;
			
				my $out = `$cmd`;
				my @res = grep(/^(NSS|Max Chi|PHI)/, split /[\n\r]/, $out);
				
				if(@res){			# tests found in output
					foreach (@res) {
						s/://;
						s/ (C|\()/_$1/;					
						my @l = split / +/;
				
						#push @res_all, [$lcb_cnt, $i+1, $i+$window, @l[0..1]]; 
						$res_all{$infile}{$l[0]} = $l[1]; 
						}
					}
				else{
						#push @res_all, [$lcb_cnt, $i+1, $i+$window, "PHI_(Normal)", 'NA'];
					$res_all{$infile}{"PHI_(Normal)"} = 'NA';
					if($all_tests){
						foreach my $q ("NSS", "Max_Chi^2"){
								#push @res_all, [$lcb_cnt, $i+1, $i+$window, $q, 'NA'];
							$res_all{$infile}{$q} = 'NA';
							}
						}
					if($phi_perm){
							#push @res_all, [$lcb_cnt, $i+1, $i+$window, "PHI_(Permutation)", 'NA'];
						$res_all{$infile}{"PHI_(Permutation)"} = 'NA';
						}
					}
				}
		};
	}
waitall;
#print Dumper %res_all;

# writing table #
foreach my $infile (sort keys %res_all){
	foreach my $test (sort keys %{$res_all{$infile}}){
		print join("\t", $infile, "NA", "NA", $test, $res_all{$infile}{$test}), "\n";
		}
	}


__END__

=pod

=head1 NAME

Phi_gene_batch.pl -- Run Phi (PhiPack) on 1 or more nucleotide alignments (fasta format)

=head1 SYNOPSIS

Phi_gene_batch.pl [options] *fasta > phi_res.txt

=head2 options

=over

=item -o  <bool>

All Phi tests (Phi, NSS, Max Chi^2). [TRUE]

-item -w  <int>

Phi window size ('-w' in Phi). [100]

=item -p  <int>

Phi permutation number ('-p' in Phi). [1000]

=item -f  <int>

Number of alignments to process in parallel. [1]

=item -v  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc Phi_gene_batch.pl

=head1 DESCRIPTION

Run Phi from PhiPack on >= 1 alignment in fasta format.

"NA" for pvalue indicates that there was not
enough phylogenetic information in the window.

Max Chi^2 and NSS can be run by using '-o'.

"start" and "end" in the output are just placeholders
needed for downstream analysis.

=head1 EXAMPLES

=head2 Basic usage

Phi_gene_batch.pl *fasta > phi_out.txt

=head2 Using 24 CPUs

Phi_gene_batch.pl -f 24 *fasta > phi_out.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

