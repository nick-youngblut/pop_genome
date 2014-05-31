#!/usr/bin/env perl

=pod

=head1 NAME

circos_featDensity.pl -- sliding window of feature density across the genome 

=head1 SYNOPSIS

circos_featDensity.pl [flags] < circos_xy.txt > featDensity_xy.txt

=head2 Required Flags

=over

=item -size  <int>

Genome size (bp).

=back

=head2 Optional Flags

=over

=item -window  <int>

Window size (bp). [10000]

=item -jump  <int>

Jump size (bp). [1000]

=item -cutoff  <float>

Value cutoff to include feature ('-compare' for sign: '>' '<'). [0]

=item -compare  <char>

Sign used for -cutoff comparisons ('>', '>=', '<', '<='). ['>=']

=back

=head2 For more information:

perldoc circos_featDensity.pl

=head1 DESCRIPTION

Get number of features within each window.

Input is a table in Circos xy format (scaffold\tstart\tstop\tvalue).

Output is in Circos xy table format.

All features (rows in table) < '-cutoff' will be skipping in the
density calculation.

=head1 EXAMPLES

=head2 Density of 4 Mb genome

circos_featDensity.pl < feat_xy.txt -s 4000000 > featDensity.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/annotation/

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
use Set::IntervalTree;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $genome_size);
my $window = 10000;
my $jump = 1000;
my $cutoff = 0;
my $compare = ">=";
GetOptions(
	   "window=i" => \$window,
	   "jump=i" => \$jump,
	   "size=i" => \$genome_size, 		# complete size of genome; must be closed 
	   "cutoff=f" => \$cutoff,
	   "compare=s" => \$compare, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults ###
die " ERROR: provide the size of the genome (bp)!\n"
	unless defined $genome_size;
check_compare($compare);

### Main ###
my $itrees_r = table2itree($cutoff, $compare);
featDensity($itrees_r, $window, $jump, $genome_size);

### subroutines ###
sub check_compare{
	my ($compare) = @_;
	die "ERROR: '-compare' must be '>', '>=', '<', or '<='\n"
		unless $compare eq '>'
		or $compare eq '>='
		or $compare eq '<'
		or $compare eq '<=';
	}	

sub featDensity{
# feature density #
	my ($itrees_r, $window, $jump, $genome_size) = @_;
	
	foreach my $scaf (keys %$itrees_r){
		for(my $i=0; $i<=$genome_size; $i+=$jump){
			my $res = $itrees_r->{$scaf}->fetch($i, $i+$window);
			print join("\t", $scaf, $i, $i+$window, scalar @$res), "\n";
			}
		}
	}

sub table2itree{
# parsing gff and placing features into an interval tree
	my ($cutoff, $compare) = @_;
	
	my %itrees;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		# sanity checks #
		die "ERROR: table should be >= 4 columns (scaffold,start,end,value)!\n"
			unless scalar @l >= 4;
		die "ERROR: 2nd & 3rd columns are not integers (should be start-end values)\n"
			unless $l[1] =~ /^\d+$/ && $l[2] =~ /^\d+$/;
		
		# applying cutoffs #
		if($compare eq '>'){
			next unless $l[3] > $cutoff;		# value > cutoff
			}
		elsif($compare eq '>='){
			next unless $l[3] >= $cutoff;
			}
		elsif($compare eq '<'){
			next unless $l[3] < $cutoff;
			}
		elsif($compare eq '<='){
			next unless $l[3] <= $cutoff;
			}
		
		# adding to interval tree #
		$itrees{$l[0]} = Set::IntervalTree->new unless exists $itrees{$l[0]};
		
		$itrees{$l[0]}->insert(1, @l[1..2]);
		}
	
	return \%itrees;
	}


