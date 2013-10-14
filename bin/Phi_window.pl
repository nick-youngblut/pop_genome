#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempfile tempdir /;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $maf_in);
my $all_tests;
my $phi_window = 100;
my $phi_perm = 1000;
my $window = 1000;
my $step = 100;
my $mult = 0;
my $fork = 0;
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
	my $pm = new Parallel::ForkManager($fork);
	
	for (my $i=0; $i<=($lcb_len-1); $i+=$step){
		$pm->start and next;
		
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
		
		my $out = `$cmd`;
		my @res = grep(/^(NSS|Max Chi|PHI)/, split /[\n\r]/, $out);
		
		if(@res){
			foreach (@res){
				s/://;
				s/ (C|\()/_$1/;
				my @l = split / +/;
			
				print join("\t", $lcb_cnt, $i+1, $i+$window, @l[0..1]), "\n"; 
				}
			}
		else{
			print join("\t", $lcb_cnt, $i+1, $i+$window, "PHI_(Normal)", 'NA'), "\n"; 
			if($all_tests){
				foreach my $q ("NSS", "Max_Chi^2"){
					print join("\t", $lcb_cnt, $i+1, $i+$window, $q, 'NA'), "\n"; 
					}
				}
			if($phi_perm){
				print join("\t", $lcb_cnt, $i+1, $i+$window, "PHI_(Permutation)", 'NA'), "\n"; 			
				}
			}
			
		$pm->finish; # do the exit in the child process
		}
	$pm->wait_all_children;
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

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

