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

Bootstrap cutoff for defining supported quartets (> -boot). [70]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc quartetDecomp_byPop.pl

=head1 DESCRIPTION

summarize output from Quartet Decomposition Server analysis

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

my ($verbose_b, $pop_in, $matrix_in, $quartets_in);
my $boot_cut = 70;		# >= bootstrap support for quartet
GetOptions(
	"population=s" => \$pop_in,
	"matrix=s" => \$matrix_in,
	"quartets=s" => \$quartets_in,
	"bootstrap=i" => \$boot_cut,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide a population file!\n" unless $pop_in;
die "ERROR: provide a matrix file!\n" unless $matrix_in;
die "ERROR: provide a quartet file!\n" unless $quartets_in;
map{ "ERROR: cannot find: '$_'\n" unless -e } ($pop_in, $matrix_in, $quartets_in);


#--- MAIN ---#
my $pops_r = load_pop($pop_in);
my $matrix_r = load_matrix($matrix_in, $boot_cut);
load_quartets($quartets_in, $pops_r, $matrix_r);

#--- Subroutines ---#
sub load_quartets{
# loading quartet index
## if quartet supported in any gene tree (matrix):
#		determine whether incongruent w/ population 
#			both quartet nodes contain both populations
	my ($quartets_in, $pops_r, $matrix_r) = @_;
	
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

			if( scalar @n1 + scalar @n2 == 4 &&									# different at each node
				scalar @n3 + scalar @n4 + scalar @n5 + scalar @n6 == 6 ){		# same across 2 pairings
				print join("\t", 
						"q\_$l[0]",								# quartet index
						$l[1],
						scalar @{$matrix_r->{$l[0]}},			# N_gene_families
						join(",", @{$matrix_r->{$l[0]}}) 		# gene_family_list
						), "\n"; 		
						
				}
				
				
			#if($l[0] == 8157){ print Dumper @n1, @n2; exit;} 
			}
		
		}
	close IN;

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
		else{
			for my $i (1..$#l){
				push @{$matrix{ $q_line{$i} }}, $l[0]
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



