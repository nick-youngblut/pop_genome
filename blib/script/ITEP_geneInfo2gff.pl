#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my @atts = ("14,note", "10,product", "1,name");
GetOptions(
		"attribute=s{,}" => \@atts,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### Main
my $atts_r = parse_atts(\@atts); 
geneInfo2gff($atts_r);


### Subroutines
sub geneInfo2gff{
#5,1,'gene',6,7,".",8,9,note="14"
	my ($atts_r) = @_;
	
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		# attributes #
		my @atts;
		foreach my $col (keys %$atts_r){
			die " ERROR! cannot find specified attribute column ", $col + 1, "\n"
				unless defined $l[$col];
			push @atts, join("=", $atts_r->{$col}, "\"$l[$col]\"");
			}
		
		# editing values #
		$l[1] =~ s/ /_/g;
		
		# prining line #
		print join("\t", $l[4], $l[1], 'CDS', @l[5..6], ".", @l[7..8], join(";", @atts)), "\n";
		}
	}

sub parse_atts{
	my ($atts_r) = @_;
	
	my %atts;
	foreach (@{$atts_r}){
		my @l = split /,/;
		#print Dumper @l;
		die " ERROR: attributes must be comma delimited pairs (column_num, attribute_label)!\n"
			unless scalar @l == 2 && $l[0] =~ /^\d+$/;
		$atts{$l[0]-1} = $l[1];
		}
	
	return \%atts;
	}

__END__

=pod

=head1 NAME

ITEP_geneInfo2gff.pl -- ITEP gene info table to a legit gff file format

=head1 SYNOPSIS

ITEP_geneInfo2gff.pl [options] < ITEPgeneInfo.txt > ITEPgeneInfo.gff

=head2 Options

=over

=item -attribute  <char>

GeneInfo table columns to add to the attributes gff column. (comma-separated pairs; ex: -a 14,note 2,taxon_name).

=item -verbose  <bool>


Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfo2gff.pl

=head1 DESCRIPTION

Convert the pseudo-gff format from ITEP to a legit gff file format.

=head1 EXAMPLES

=head2 Basic usage:

ITEP_geneInfo2gff.pl < ITEPgeneInfo.txt > ITEP.gff

=head2 ClusterID as note in attribute column & taxon_name 

ITEP_geneInfo2gff.pl -a 14,note 2,taxon_name < ITEPgeneInfo.txt > ITEP.gff

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITOL_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

