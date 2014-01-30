#!/usr/bin/env perl

=pod

=head1 NAME

quartetDecomp_distCut.pl -- filter out quartets that do not have differentiating SNPs

=head1 SYNOPSIS

quartetDecomp_distCut.pl [flags] < quartet_summary.txt > quartet_summary_filtered.txt

=head2 Required flags

=over

=item -align  <char>

Alignment file (2-column; tab-delimited; file\\tID). 

=back

=head2 Optional flags

=over

=item -SNP  <int>

Number of SNPs that must differentiate taxa on either node of quartet. [1]

=item -regex  <char>

Regex to alter taxa names in alignment files (2 arguments required; eg. -r '.peg.+' ''). 

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc quartetDecomp_distCut.pl

=head1 DESCRIPTION

Some quartets that passed the bootstrap support cutoff
may actually have no differentiating SNPs, so the 
supported topology was just by chance. This script will
remove any such quartets (for each tree/alignment) 
that do not meet the specified SNP requirement.

The ouput will be by ID. Columns include: 'ID'\t'quartetID'\t'quartet_tree'

=head1 EXAMPLES

=head2 Basic usage:

quartetDecomp_distCut.pl -a align_list.txt < quartet_summary.txt > quartet_summary_filtered.txt

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

my ($verbose_b, $align_list, @regex);
my $snp_cut = 1;
GetOptions(
	"snp=i" => \$snp_cut, 
	"align=s" => \$align_list,
	"regex=s{2,2}" => \@regex,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide an alignment list file!\n"
	unless defined $align_list;
die "ERROR: cannot find $align_list!\n"
	unless -e $align_list;
$regex[0] = qr/$regex[0]/ if defined $regex[0];	

#--- MAIN ---#
my $align_list_r = load_align_list($align_list);
my $tbl_r = load_quartet_by_cluster();
check_dists($tbl_r, $align_list_r, \@regex);
write_tbl_byCluster($tbl_r);


#--- Subroutines ---#
sub write_tbl_byCluster{
# writing out editted table: by Cluster #
	my ($tbl_r) = @_;

	my $clust_cnt = 0;
	foreach my $clust (sort keys %$tbl_r){
		foreach my $quartet (sort keys %{$tbl_r->{$clust}}){
			print join("\t", $quartet, 
					$tbl_r->{$clust}{$quartet}{'status'},
					$tbl_r->{$clust}{$quartet}{'tree'},
					$clust), "\n";
			}
		$clust_cnt++ if scalar keys %{$tbl_r->{$clust}};
		}
		
	print STDERR "Number of remaining clusterIDs after filtering: $clust_cnt\n";
	}

sub check_dists{
	my ($tbl_r, $align_list_r, $regex_r) = @_;
	
	my $del_cnt = 0;
	foreach my $clust (keys %$tbl_r){		# all clusters w/ supported quartets
		# status #
		print STDERR "Processing cluster: $clust\n" unless $verbose_b;
		
		# loading fasta #
		die "ERROR: cannot find '$clust' in alignment list IDs\n"
			unless exists $align_list_r->{$clust};
		my $fasta_r = load_fasta($align_list_r->{$clust}, $regex_r); 
		
		foreach my $quartet (keys %{$tbl_r->{$clust}}){
			my $n_var = get_pairwise_snp($fasta_r, $tbl_r->{$clust}{$quartet}{'taxa'}, $snp_cut);
			
			# should have n_var of 4 if all vary #
			if($n_var == 4){ 
				next;
				}
			elsif($n_var < 4){ 
				warn "Only $n_var of the 4 across-quartet comparisons in ($clust -> $quartet) have >=$snp_cut SNPs. Deleting the gene_fam->quartet!\n"
					unless $verbose_b;
				delete $tbl_r->{$clust}{$quartet};
				$del_cnt++;
				}
			else{ die "LOGIC ERROR: $!\n"; }
			}
		}
		
	# status #
	print STDERR "Number of quartets deleted from any clusterID: $del_cnt\n";
		#print Dumper %$tbl_r; exit;
	}
	
sub get_pairwise_snp{
# getting number of SNPs when comparing each taxon in each side of quartet #
	my ($fasta_r, $taxa_r, $snp_cut) = @_;
	
	# getting N-snps on pairwise comparisons of snps on either side of quartet #
	## all pairwise comparisons of members on each side of quartet ##
	my $n_var = 0; 				# must have all 4 comparisons w/ sequence variation
	for my $i (0..1){			# quartet side 1
		for my $ii (0..1){		# quartet side 2
			my $n_snp = get_snp(
					$$taxa_r[0][$i],
					$$taxa_r[1][$ii],
					$fasta_r->{$$taxa_r[0][$i]}, 
					$fasta_r->{$$taxa_r[1][$ii]}
					);
			$n_var++ if $n_snp >= $snp_cut;
			}
		}
	return $n_var;
	}

sub get_snp{
	my ($tax1, $tax2, $seq1, $seq2) = @_;
	die "ERROR: $tax1 sequence not found!\n"
		unless defined $seq1;
	die "ERROR: $tax2 sequence not found!\n"
		unless defined $seq2;
	die "ERROR: $tax1 & $tax2 sequences are not the same length!\n"
		unless scalar @$seq1 == scalar @$seq2;
	
	my $snp_cnt = 0;
	for my $i (0..($#$seq1)){
		if( $$seq1[$i] ne $$seq2[$i] 
			&& $$seq1[$i] ne '-'
			&& $$seq2[$i] ne '-'){
			$snp_cnt++;
			}	
		}
		
	return $snp_cnt;		# N_snps
	}
	
sub load_fasta{
# loading fasta file as a hash #
	my ($fasta_in, $regex_r) = @_;
	
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if  /^\s*$/;	
 		
 		if(/>.+/){
 			($tmpkey = $_) =~ s/^>//;	# changing key
 			$tmpkey =~ s/$$regex_r[0]/$$regex_r[1]/ if defined $regex_r;
 			}
 		else{
 			push @{$fasta{$tmpkey}}, split //, $_;
 			}
		}
	close IN;

		#rint Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta

sub load_quartet_by_cluster{
# loading quartet list (output of quartetDecomp_byPop) by cluster #
	
	my %tbl;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die "ERROR: quartet summary should have 4 columns!\n"
			unless scalar @l == 4;
		my @ll = split /,/, $l[3];
		
		# quartet tree #
		my @q = split /,/, $l[2];
		map{ s/[();]+//g } @q;
		
		# loading hash #
		## cluster_ID=>quartetID=>cat=>value			
		$tbl{$l[3]}{$l[0]}{'tree'} = $l[2];
		$tbl{$l[3]}{$l[0]}{'taxa'} = [ [@q[0..1]], [@q[2..3]] ];
		$tbl{$l[3]}{$l[0]}{'status'} = $l[1];
		}

	# status #
	print STDERR "Number of clusterIDs in quartet summary file: ",
		scalar keys %tbl, "\n";
		#print Dumper %tbl; exit;
	return \%tbl;
	}

sub load_align_list{
	my ($align_list) = @_;

	open IN, $align_list or die $!;	
	my %align_list;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		die "ERROR: alignment list is not in 2 column format: file\\tID\n"
			unless scalar @l == 2;
		warn "WARNING: '$l[1]' is duplicated in alignment list! Just using last instance\n"
			if exists $align_list{$l[1]};
		$align_list{$l[1]} = $l[0];			# ID => file
		}
	close IN;
	
		#print Dumper %align_list; exit;
	return \%align_list;
	}
	
