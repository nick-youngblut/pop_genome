#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/min/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $write_NA, $make_neg);
GetOptions(
	   "NA" => \$write_NA, 
	   "strand" => \$make_neg,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $info_r = load_gene_info_table();
find_distance($info_r);

### Subroutines
sub find_distance{
# finding distance between genes on same scaffold #
	my ($info_r) = @_;
	foreach my $taxon_name (sort keys %$info_r){
		
		# figs on same scaffold in same taxon #
		my @fig = keys %{$info_r->{$taxon_name}};
		for my $i (0..$#fig){
			for my $ii (0..$#fig){
				next if $ii <= $i;
				
				# same scaffold? #
				if(${$info_r->{$taxon_name}{$fig[$i]}}[4] ne 
					${$info_r->{$taxon_name}{$fig[$ii]}}[4]){
					
					print join("\t", "NA", 
						@{$info_r->{$taxon_name}{$fig[$i]}},
						@{$info_r->{$taxon_name}{$fig[$ii]}}), "\n"
						if $write_NA;
					}
				
				# start - end #
				my $start1 = ${$info_r->{$taxon_name}{$fig[$i]}}[5];
				my $end1 = ${$info_r->{$taxon_name}{$fig[$i]}}[6];
				my $start2 = ${$info_r->{$taxon_name}{$fig[$ii]}}[5];
				my $end2 = ${$info_r->{$taxon_name}{$fig[$ii]}}[6];
				
				# finding min distance #
				my @dist;
				push(@dist, abs($start1 - $start2));
				push(@dist, abs($end1 - $start2));
				push(@dist, abs($start1 - $end2));					
				push(@dist, abs($end1 - $end2));					
				
				my $min_dist = min @dist;
				
				# making negative if different strand #
				$min_dist = -$min_dist 
					if ${$info_r->{$taxon_name}{$fig[$i]}}[8] 
						ne ${$info_r->{$taxon_name}{$fig[$ii]}}[8]
					&& $make_neg;
				
				# writing out lines #
				print join("\t", $min_dist, 
					@{$info_r->{$taxon_name}{$fig[$i]}},
					@{$info_r->{$taxon_name}{$fig[$ii]}}), "\n";
				}
			}
		}	
	}

sub load_gene_info_table{
# loading gene info table #
## taxon_name => fig => row
	my %info;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		$info{$line[1]}{$line[0]} = \@line;
		}

		#print Dumper %info; exit;
	return \%info;
	}


__END__

=pod

=head1 NAME

ITEP_geneInfo2scaffoldDist.pl -- getting distance between genes on same scaffold

=head1 SYNOPSIS

ITEP_geneInfo2scaffoldDist.pl [options] < geneInfoTable.txt > distanceInfoTable.txt

=head2 options

=over

=item -NA  <bool>

Write gene comparisons if they were not on the same scaffold ('NA' for distance)? [FALSE]

=item -strand  <bool>

Make distances negative if genes are on different strands? [FALSE]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfo2scaffoldDist.pl

=head1 DESCRIPTION

Determine minimum distance (bp) between genes found
in the same taxon on the same scaffold.

=head2 Output

The 1st column is minimum distance (bp)
The following columns are the original rows 
for both genes appended together.

=head1 EXAMPLES

=head2 Basic usage:

db_getGenesWithAnnotation.py "methyl coenzyme M reductase" | db_getGeneInformation.py | 
ITEP_geneInfo2scaffoldDist.pl

=head2 Just taxon, distance, geneIDs, and annotations

db_getGenesWithAnnotation.py "methyl coenzyme M reductase" | db_getGeneInformation.py | 
ITEP_geneInfo2scaffoldDist.pl | cut -f 1-3,11,23 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

