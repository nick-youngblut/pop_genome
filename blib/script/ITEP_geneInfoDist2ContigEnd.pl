#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use List::Util qw/min/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $scaf_b);
my $org_col = 3;
GetOptions(
	   "organism=i" => \$org_col,
	   "scaffold" => \$scaf_b, 		# don't split scaffolds into contigs? [fALSE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$org_col--;


### Main
my $gi_r = load_geneInfo($org_col);
my $org_list_file = write_org_list($gi_r);
my $contig_r = get_contigs($gi_r, $org_list_file);
get_min_dist2scaffold_end($gi_r, $contig_r);

### subroutines
sub get_min_dist2scaffold_end{
	my ($gi_r, $contig_r) = @_;
	
	foreach my $orgID (keys %$gi_r){
		foreach my $line (@{$gi_r->{$orgID}}){
			# getting min distanct between any gap or scaffold start-end #
			my $min = min (@{$line}[5..6], $contig_r->{$$line[4]});
			print join("\t", @$line, $min), "\n";
			}
		}
	}

sub get_contigs{
# calling db_getContigs.py #
	my ($gi_r, $org_list_file) = @_;

	die " ERROR: cannot find $org_list_file!\n"
		unless -e $org_list_file;
	open PIPE, "cat $org_list_file | db_getContigs.py | db_getContigs.py | db_getContigSeqs.py | "
		or die $!;

	my %contig;
	# loading locations of gaps and ends #
	while(<PIPE>){
		chomp;
		my @l = split /\t/; 

		# getting locations of gaps #
		unless($scaf_b){
			my $i = 0;
			while($l[1] =~ /N+/gi){
				$i++;
				#$contig{$l[0]}{"gap$i"}{"start"} = $-[0] + 1;
				#$contig{$l[0]}{"gap$i"}{"end"} = $+[0] + 1;
				push @{$contig{$l[0]}}, $-[0] + 1;
				push @{$contig{$l[0]}}, $+[0] + 1;
				}
			}
		# scaffold start-end #
		#$contig{$l[0]}{"scaf"}{"start"} = 1;
		#$contig{$l[0]}{"scaf"}{"stop"} = length($l[1]);
		push @{$contig{$l[0]}}, 1;
		push @{$contig{$l[0]}}, length($l[1]);
		}
	close PIPE;
	unlink $org_list_file or die $!;
		
		#print Dumper %contig; exit;
	return \%contig;		# start-stop 1-indexed
	}
	
sub write_org_list{
	my ($gi_r) =@_;
	
	my $outfile = "tmp_org_list.txt";
	open OUT, ">$outfile" or die $!;
	foreach my $orgID (keys %$gi_r){
		print OUT $orgID, "\n";
		}
	close OUT;
	
	return $outfile;
	}

sub load_geneInfo{
# loading geneInfo; getting unique organism ID list [3rd column] #
	my ($org_col) = @_;
	my %gi;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die " ERROR organism column not found!\n" 
			unless defined $l[$org_col];
		
		push @{$gi{$l[$org_col]}}, \@l;
		}
		
		#print Dumper %gi; exit;
	return \%gi;
	}



__END__

=pod

=head1 NAME

ITEP_geneInfoDist2ContigEnd.pl -- determine the minimum distince between a gene and a contig end

=head1 SYNOPSIS

ITEP_geneInfoDist2ContigEnd.pl [options] < geneInfo.txt > geneInfo_minDist.txt

=head2 Flags

=over

=item -org  <int>

Column in the gene info table containing the organism ID. 

=item -scaf  <bool>

Min distance to end of the scaffold instead of end of a contig? [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfoDist2ContigEnd.pl

=head1 DESCRIPTION

Determine whether a gene is at the end of a contig/scaffold.
The input is an ITEP geneInfo table.
The minimum distance (bp) between a gene start or end
to a contig/scaffold start end is appended to the end 
of the output geneInfo table. 

=head1 EXAMPLES

=head2 Basic usage:

ITEP_geneInfoDist2ContigEnd.pl < geneInfo.txt > geneInfo_minDist.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

