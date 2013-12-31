#!/usr/bin/env perl

=pod

=head1 NAME

quartetDecomp_byPop.pl -- summarize output from Quartet Decomposition Server analysis

=head1 SYNOPSIS

quartetDecomp_byPop.pl [flags] > quartet_summary.txt

=head2 Required flags

=over

=item -pop  <char>

Population file (2-column; tab-delimited; taxon\tpopulation)

=item -matrix  <char>

matrix.txt file from Quartet Decomposition Server

=item -quartet  <char>

quartets.txt file from Quartet Decomposition Server

=back

=head2 Optional flags

=over

=item -boot  <int>

Bootstrap cutoff for defining supported quartets (> -boot). [80]

=item -gene  <char>

genefamliy.out file from Quartet Decomposition Server. (used to label gene families)

-item -regex  <char>

Regex to alter file names in genefamliy.out (2 arguments required; eg. -r '.+maxbit|\.fasta_ML' ''). 

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc quartetDecomp_byPop.pl

=head1 DESCRIPTION

summarize output from Quartet Decomposition Server analysis

=head2 output

=over

=item * quartet ID

=item * 'cogruent' or 'congruent' depending on quartet topology

=item * quartet tree

=item * gene tree containing quartet (modified by -regex if provided)

=back

=head1 EXAMPLES

=head2 Basic usage:

quartetDecomp_byPop.pl [flags] > quartet_summary.txt

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

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $pop_in, $matrix_in, $quartets_in, $gene_fam_in, @regex);
my $boot_cut = 80;		# >= bootstrap support for quartet
GetOptions(
	"population=s" => \$pop_in,
	"matrix=s" => \$matrix_in,
	"quartets=s" => \$quartets_in,
	"bootstrap=i" => \$boot_cut,
	"gene=s" => \$gene_fam_in,
	"regex=s{2,2}" => \@regex,	
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide a population file!\n" unless $pop_in;
die "ERROR: provide a matrix file!\n" unless $matrix_in;
die "ERROR: provide a quartet file!\n" unless $quartets_in;
map{ "ERROR: cannot find: '$_'\n" unless -e } ($pop_in, $matrix_in, $quartets_in);
$regex[0] = qr/$regex[0]/ if defined $regex[0];

#--- MAIN ---#
my $pops_r = load_pop($pop_in);
my $gene_fam_r = load_gene_fam($gene_fam_in, \@regex) if $gene_fam_in;
my $matrix_r = load_matrix($matrix_in, $boot_cut);
load_quartets($quartets_in, $pops_r, $matrix_r, $gene_fam_r);

#--- Subroutines ---#
sub load_quartets{
# loading quartet index
## if quartet supported in any gene tree (matrix):
#		determine whether incongruent w/ population 
#			both quartet nodes contain both populations
	my ($quartets_in, $pops_r, $matrix_r, $gene_fam_r) = @_;
	
	open IN, $quartets_in or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die "ERROR: quartet file should have 2 columns: quartet_index\\tquartet tree!\n"
			unless scalar @l == 2;
		
		if(exists $matrix_r->{$l[0]}){		# quartet supported in >= 1 gene tree
			# parsing quartet tree #
			(my $tmp = $l[1]) =~ s/[();]+//g;
			my @taxa = split /,/, $tmp;
			die "ERROR: '$_' does not have 4 taxa!\n"
				unless scalar @taxa == 4;
		
			# determining whether mixed population on each quartet node #
			map{ die "ERROR: cannot find $_ in population file!\n" unless exists
					$pops_r->{$_} } @taxa;
					
			my @n1 = unique( $pops_r->{$taxa[0]}, $pops_r->{$taxa[1]} );
			my @n2 = unique( $pops_r->{$taxa[2]}, $pops_r->{$taxa[3]} );
			my @n3 = unique( $pops_r->{$taxa[0]}, $pops_r->{$taxa[2]} );
			my @n4 = unique( $pops_r->{$taxa[1]}, $pops_r->{$taxa[3]} );
			my @n5 = unique( $pops_r->{$taxa[0]}, $pops_r->{$taxa[3]} );
			my @n6 = unique( $pops_r->{$taxa[1]}, $pops_r->{$taxa[2]} );

			# converting gene family to user-defined name (if provided) #
			if($gene_fam_r){
				foreach my $fam (@{$matrix_r->{$l[0]}}){
					$fam =~ s/^fam_//;
					die "ERROR: cannot find '$fam' in gene family file!\n"
						unless exists $gene_fam_r->{$fam};
					$fam = $gene_fam_r->{$fam};
					}
				}
			
			# writing out line #
			my $status; 
			if( scalar @n1 + scalar @n2 == 4 &&									# different at each node
				scalar @n3 + scalar @n4 + scalar @n5 + scalar @n6 == 6 ){		# same across 2 pairings
				$status = 'incongruent';
				}
			else{ $status = 'congruent'; }
			
			foreach my $gene ( @{$matrix_r->{$l[0]}} ){
				print join("\t", 
					"q_$l[0]",
					$status,						
					$l[1],
					$gene
					), "\n"; 		
				}
			}
		}
	close IN;

	}

sub load_gene_fam{
	my ($gene_fam_in, $regex_r) = @_;
	
	open IN, $gene_fam_in or die $!;
	my %gene_fam;
	while(<IN>){
		chomp;
		my @l = split /\t/;
		die "ERROR: gene family order file should be 3 columns: order\\tfile_name\\tN_trees\n"
			unless scalar @l == 3 && $l[0] =~ /^\d+$/;
		
		$l[1] =~ s/$$regex_r[0]/$$regex_r[1]/g if defined $$regex_r[0];
			
		die "ERROR: duplicate ordering values in gene family order list!\n"
			if exists $gene_fam{$l[0]};
		$gene_fam{$l[0]} = $l[1];
		}
	close IN;
	
		#print Dumper %gene_fam; exit;
	return \%gene_fam;
	}
	
sub unique{
	my %u;
	map{ $u{$_} = 1 } @_;
	return keys %u;
	}
	
sub load_matrix{
# loading matrix & which trees (gene_families) support the quartet #
	my ($matrix_in, $boot_cut) = @_;
	
	open IN, $matrix_in or die $!;
	my %q_line;
	my %matrix;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		if($. == 1){
			for my $i (1..$#l){
				$l[$i] =~ s/^q_//;			# just quartet index 
				$q_line{$i} = $l[$i];
				}
			}
		else{			# applying bootstrap cutoff 
			for my $i (1..$#l){
				push @{$matrix{ $q_line{$i} }}, $l[0]			# all gene families w/ supported quartet
					if $l[$i] > $boot_cut;
				}
			}
		}
	close IN;
		#print Dumper %matrix; exit;	
	return \%matrix;
	}

sub load_pop{
# loading population file designating taxon=>population #
	my ($pop_in) = @_;
	
	open IN, $pop_in or die $!;
	my %pop;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die "ERROR: population file should have 2 columns: (taxon_name\\tpopulation)!\n"
			unless scalar @l == 2;
		warn "WARNING: $l[0] duplicated in population file! Only last instance used!\n"
			if exists $pop{$l[0]};
		$pop{$l[0]} = $l[1];
		}
	close IN;

	return \%pop;
	}



