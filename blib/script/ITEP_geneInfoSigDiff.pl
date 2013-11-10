#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
 use Statistics::Test::WilcoxonRankSum;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $pop_in);
GetOptions(
	   "pop=s" => \$pop_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a population file (name\tpopulation)\n"
	unless $pop_in;
die " ERROR: can't find $pop_in\n" unless -e $pop_in;

### Main
my ($pop_r, $upop_r) = load_pop($pop_in);
my $info_r = load_geneInfo($pop_r);
get_sig_diff($info_r, $upop_r, "copy_number");
get_sig_diff($info_r, $upop_r, "gene_length");

### subroutines 
sub get_sig_diff{
	my ($info_r, $upop_r, $cat) = @_;
	
	foreach my $clusterID (keys %$info_r){
		# pairwise comparisons of populations 
		for my $i (0..$#$upop_r){
			for my $ii (0..$#$upop_r){
				next if $i >= $ii;
				my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
				
				$wilcox_test->load_data( $info_r->{$clusterID}{$cat}{$$upop_r[$i]},
										$info_r->{$clusterID}{$cat}{$$upop_r[$ii]});
				print join("\t", $clusterID, $$upop_r[$i], $$upop_r[$ii], 
							scalar @{$info_r->{$clusterID}{$cat}{$$upop_r[$i]}},
							scalar @{$info_r->{$clusterID}{$cat}{$$upop_r[$ii]}},
							$cat, 
							$wilcox_test->probability()), "\n";
				}
			}
		}
	}

sub get_sig_diff_copy{
	my ($info_r, $upop_r) = @_;
	
	foreach my $clusterID (keys %$info_r){
		# pairwise comparisons of populations 
		for my $i (0..$#$upop_r){
			for my $ii (0..$#$upop_r){
				next if $i >= $ii;
				my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
				
				$wilcox_test->load_data( $info_r->{$clusterID}{"copy"}{$$upop_r[$i]},
										$info_r->{$clusterID}{"copy"}{$$upop_r[$ii]});
				print join("\t", $clusterID, $$upop_r[$i], $$upop_r[$ii], "copy_number", 
							$wilcox_test->probability()), "\n";
				}
			}
		}
	}


sub load_geneInfo{
# loading geneInfo table #
	my ($pop_r) = @_;
	
	my %info;
	my %warnings;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
	
		# sanity check #
		die " ERROR: no clusterID column in geneInfo table! You must use db_getClusterGeneInformation!\n"	
			unless scalar @l >=14 && $l[13] =~ /^\d+$/;
		unless(exists $pop_r->{$l[1]}){
			$warnings{$l[1]}=1;
			next;
			}
	
		# getting copy & length #
		$info{$l[13]}{"copy_bytax"}{$pop_r->{$l[1]}}{$l[2]}++;		# copy per taxon
		push @{$info{$l[13]}{"gene_length"}{$pop_r->{$l[1]}}}, length $l[10];		# copy per taxon
		}
	
	# getting copy #
	foreach my $clusterID (keys %info){
		foreach my $pop (keys %{$info{$clusterID}{"copy_bytax"}}){
			foreach my $taxon (keys %{$info{$clusterID}{"copy_bytax"}{$pop}}){
				push @{$info{$clusterID}{"copy_number"}{$pop}}, $info{$clusterID}{"copy_bytax"}{$pop}{$taxon};
				}
			}
		}
	
	# writing out warnings #
	foreach my $taxon (keys %warnings){
		print STDERR "WARNING: '$taxon' in 2nd column of geneInfo table not found in population table!\n"
		}
		
		#print Dumper %info; exit;
	return \%info;
	}

sub load_pop{
	my ($pop_in) = @_;
	
	open IN, $pop_in or die $!;
	
	my %pop;
	my %upop;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die " ERROR: pop table must have 2 columns (taxon	pop)\n"
			unless scalar @l >= 2;
			
		$pop{$l[0]} = $l[1];
		$upop{$l[1]} = 1;
		}
	close IN;
	
	die " ERROR: >=2 populations required!\n"
		unless scalar keys %upop >= 2;
	
	return \%pop, [sort keys %upop];
	}



__END__

=pod

=head1 NAME

ITEP_geneInfoSigDiff.pl -- Check for signif diffs in gene copy number of length between populations

=head1 SYNOPSIS

ITEP_geneInfoSigDiff.pl [options] < ITEPgeneInfo.txt > sig-diff_summary.txt

=head2 Options

=over

=item -pop  <char>

Population table (taxonID\tpopulationID). [required]

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfoSigDiff.pl

=head1 DESCRIPTION

Mann-Whitney U test on gene copy number and gene length
for each cluster. 

You need >= 2 populations. 

=head2 Output

=over

=item clusterID

=item pop1

=item pop2

=item N_pop1

=item N_pop2

=item category

=item P_value

=back

=head1 EXAMPLES

=head2 Basic usage:

ITEP_geneInfoSigDiff.pl < ITEPgeneInfo.txt -p pop.txt > sig-diff_summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

