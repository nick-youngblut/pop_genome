#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempfile tempdir /;
#use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $maf_in);
my $all_tests;
my $phi_window = 100;
my $phi_perm = 1000;
my $window = 1000;
my $step = 100;
my $mult = 0;
my $fork = 1;
GetOptions(
	   "maf=s" => \$maf_in,
	   "n=i" => \$window,
	   "step=i" => \$step,
	   "o" => \$all_tests,			# TRUE
	   "mult=i" => \$mult,
	   "w" => \$phi_window,
	   "p" => \$phi_perm,
	   "fork=i" => \$fork, 		
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a maf file!\n"
	unless $maf_in;
die " ERROR: cannot find $maf_in!\n"
	unless -e $maf_in;

### MAIN
# tempdir #
my $tmpdir = File::Temp->newdir();
my $dirname = $tmpdir->dirname;


# phi by LCB #
phi_by_lcb($maf_in, $dirname);


### Subroutines
sub phi_by_lcb{
# run phi along a sliding window and keep track of the LCB; postions by LCB #
	my ($maf_in, $dirname) = @_;
	
	# output header #
	print join("\t", qw/LCB start end test pvalue/), "\n";
	
	# loading maf #
	open IN, $maf_in or die $!;
	my $skip= 0;
	my $lcb_cnt = 0;
	my %lcb_seq;
	my $lcb_len;
	while(<IN>){
		chomp;
		if(/^a /){
			$lcb_cnt++;
			my @l = split / /;			
			(my $m = $l[3]) =~ s/mult=//;
			$skip = 1 if $m < $mult;
			
			# status #
			if($skip){
				print STDERR " Skipping LCB$lcb_cnt, number of taxa in LCB is < '-mult'\n";
				}
			else{ 
				print STDERR " Processing LCB$lcb_cnt\n" unless $verbose; 
				}
			}
		elsif(/^s /){
			next if $skip;		# skipping if necessary
			my @l = split /[\t ]+/;
			(my $taxon = $l[1]) =~ s/\..+//;
			$lcb_seq{$taxon} = $l[6];
				
			$lcb_len = length $l[6];
			
			}
		elsif(/^\s*$/ ){		# end of LCB; calling phi; resetting;  			
			$skip = 0;
			next unless %lcb_seq;			# skipping unless sequences found
			
			# calling phi #
			call_phi(\%lcb_seq, $lcb_cnt, $lcb_len, $dirname);
			
			%lcb_seq = ();
			}
		}
	}


sub call_phi{
# for sliding window, writing of
	my ($seqs_r, $lcb_cnt, $lcb_len, $dirname) = @_;
	
	# forking #
	use Forks::Super;
	
	my %res_all;
	for (my $i=0; $i<=($lcb_len-1); $i+=$step){
	#for (my $i=0; $i<=10000; $i+=$step){
		#$pm->start and next;
		my $job = fork {
			max_proc => $fork,
			share => [ \%res_all ],
			sub => sub {
				#push @res, "here";
				# writing out temporary file #
				my $tmp_file = "$dirname/tmp$i.fasta";
				open OUT, ">$tmp_file" or die $!;
				foreach my $taxon (keys %$seqs_r){
					print OUT join("\n", ">$taxon", substr $seqs_r->{$taxon}, $i, $window), "\n";
					}
				close OUT;
			
				# calling phi #
				my $cmd = "Phi -f $tmp_file -w $phi_window";
				$cmd .= " -p $phi_perm" if $phi_perm;
				$cmd .= " -o" if $all_tests;
				$cmd .= " 2>/dev/null" unless $verbose;
			
				my $out = `$cmd`;
				my @res = grep(/^(NSS|Max Chi|PHI)/, split /[\n\r]/, $out);
			
				if(@res){			# tests found in output
					foreach (@res) {
						s/://;
						s/ (C|\()/_$1/;					
						my @l = split / +/;
				
						#push @res_all, [$lcb_cnt, $i+1, $i+$window, @l[0..1]]; 
						$res_all{$i+1}{$l[0]} = $l[1]; 
						}
					}
				else{
						#push @res_all, [$lcb_cnt, $i+1, $i+$window, "PHI_(Normal)", 'NA'];
					$res_all{$i+1}{"PHI_(Normal)"} = 'NA';
					if($all_tests){
						foreach my $q ("NSS", "Max_Chi^2"){
								#push @res_all, [$lcb_cnt, $i+1, $i+$window, $q, 'NA'];
							$res_all{$i+1}{$q} = 'NA';
							}
						}
					if($phi_perm){
							#push @res_all, [$lcb_cnt, $i+1, $i+$window, "PHI_(Permutation)", 'NA'];
						$res_all{$i+1}{"PHI_(Permutation)"} = 'NA';
						}
					}
				}
			};
		}
	waitall;

	# writing table #
	foreach my $start (sort{$a<=>$b} keys %res_all){
		foreach my $test (sort keys %{$res_all{$start}}){
			print join("\t", $lcb_cnt, $start, $start+$window, $test, $res_all{$start}{$test}), "\n";
			}
		}
	}


__END__

=pod

=head1 NAME

Phi_window.pl -- Run Phi (PhiPack) along a sliding window

=head1 SYNOPSIS

Phi_window.pl [options] -maf file.maf > phi_out.txt

=head2 options

=over

=item -maf  <char>

maf file name.

=item -mult  <int>

LCBs with < '-mult' taxa will be skipped. [0]

=item -n  <int>

Window size (bp). [1000]

=item -s  <int>

Step size (bp). [100]

=item -o  <bool>

All Phi tests (Phi, NSS, Max Chi^2). [TRUE]

-item -w  <int>

Phi window size ('-w' in Phi). [100]

=item -p  <int>

Phi permutation number ('-p' in Phi). [1000]

=item -f  <int>

Number of parallel calls of Phi. [1]

=item -v  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc Phi_window.pl

=head1 DESCRIPTION

Run Phi from PhiPack using a sliding window.

"NA" for pvalue indicates that there was not
enough phylogenetic information in the window.

Alignment positions indexed by 1.

=head1 EXAMPLES

=head2 Basic usage

Phi_window.pl -maf file.maf > phi_out.txt

=head2 Using 24 CPUs

Phi_window.pl -f 24 -maf file.maf > phi_out.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

