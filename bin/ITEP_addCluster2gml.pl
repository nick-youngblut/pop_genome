#!/usr/bin/env perl

=pod

=head1 NAME

ITEP_addCluster2gml.pl -- adding clusterID to a gml file produced by db_makeGraphFromBlastResults.py

=head1 SYNOPSIS

db_makeGraphFromBlastResults.py | ITEP_addCluster2gml.pl [flags] > file.gml

=head2 Required flags

=over

=item -r

ITEP runID

=back

=head2 Optional flags

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_addCluster2gml.pl

=head1 DESCRIPTION

A simple script for adding ITEP cluster IDs to a gml file.

=head1 EXAMPLES

=head2 Basic usage:

ITEP_addCluster2gml.pl < in.gml -r all_I_2.0_c_0.4_m_maxbit > out.gml

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp qw/ tempdir /;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b);
my $runID;
GetOptions(
	"runID=s" => \$runID,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide an ITEP runID!\n"
	unless defined $runID;

#--- MAIN ---#
parse_gml($runID);

#--- Subroutines ---#
sub parse_gml{
# calling ITEP 1 node at a time #
	my ($runID) = @_;
	
	while(<>){
		chomp;
		print $_, "\n";
		if(/^\s+label "fig.+peg/){
			s/.+"(.+)".*/$1/;
			open PIPE, "echo \"$_\" | db_getClustersContainingGenes.py | grep $runID |" or die $!;
			while(<PIPE>){
				my @tmp = split /\t/;
				print "    cluster \"$tmp[1]\"\n"
				}
			close PIPE;
			}
		}	
	}
