#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Statistics::Descriptive;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $strand_comb);
GetOptions(
	   "strand" => \$strand_comb, 		# combine features from opposite strands? [FALSE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### Main
my %clust;
while(<>){
	chomp;
	next if /^\s*$/;
	my @l = split /\t/;
	
	# sanity check #
	die " ERROR: no clusterID column in geneInfo table! You must use db_getClusterGeneInformation!\n"	
		unless scalar @l >=14 && $l[13] =~ /^\d+$/;
	
	# cluster #
	my $clustID = $l[13];

	# getting start-stop-strand #
	my $start = $l[5];
	my $stop = $l[6];
	my $strand = $l[7];
	$strand = "+" if $strand_comb; 						# not accounting for strand (all to + strand)
	
	# flipping start & stop if needed #
	($start,$stop) = ($stop,$start) if $strand eq "-" && $strand_comb; 	# all to + strand
	

	# intializing #
	$clust{$clustID}{$strand}{"start"} = Statistics::Descriptive::Full->new()
		unless exists $clust{$clustID}{$strand}{"start"};
	$clust{$clustID}{$strand}{"stop"} = Statistics::Descriptive::Full->new()
		unless exists $clust{$clustID}{$strand}{"stop"};
		
	# loading #
	$clust{$clustID}{$strand}{"start"}->add_data($start);
	$clust{$clustID}{$strand}{"stop"}->add_data($stop);

	$clust{$clustID}{$strand}{"data"} = [@l[(0..4,8..$#l)]];
	}

# writting out lines w/ median #
foreach my $clustID (keys %clust){
	foreach my $strand (keys %{$clust{$clustID}}){
		my $median_start = int($clust{$clustID}{$strand}{"start"}->median());
		my $median_stop = int($clust{$clustID}{$strand}{"stop"}->median());
		
		print join("\t", @{$clust{$clustID}{$strand}{"data"}}[0..4],
						$median_start, $median_stop, $strand,
						@{$clust{$clustID}{$strand}{"data"}}[5..$#{$clust{$clustID}{$strand}{"data"}}]),
						"\n";
		}

	}



__END__

=pod

=head1 NAME

ITEP_geneInfoMedianCoords.pl -- ITEP gene info table to just 1 entry with median start-stops

=head1 SYNOPSIS

ITEP_geneInfoMedianCoords.pl [options] < ITEPgeneInfo.txt > ITEPgeneInfo_median.txt

=head2 Options

=over

=item -strand  <char>

Ignore differences in strand within a cluster (all entries will be + strand)? [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfoMedianCoords.pl

=head1 DESCRIPTION

Get the general median start-stops of
genes in
a

=head1 EXAMPLES

=head2 Basic usage:

ITEP_geneInfoMedianCoords.pl < ITEPgeneInfo.txt > ITEP.gff

=head2 ClusterID as note in attribute column & taxon_name 

ITEP_geneInfoMedianCoords.pl -a 14,note 2,taxon_name < ITEPgeneInfo.txt > ITEP.gff

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITOL_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

