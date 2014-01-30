#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my ($tree_in,$clades_in);
my $tformat = "newick";
my $brlen_cut;
my $absent_cut = 0;
my $subclade_cut = 0;
GetOptions(
		"tree=s" => \$tree_in,				# tree file
		"format=s" => \$tformat,				# tree formate
		"clades=s" => \$clades_in, 			# 2-column table: taxon\tclade
		"brlen=f" => \$brlen_cut, 			# brlen for cutoff
		"absent=f" => \$absent_cut, 		# number of taxa that can be missing a gene in cluster (if < 1, fraction) [<=]
		"subclade=f" => \$subclade_cut, 	# subclade cutoff (fraction of subclade that must have cluster) [>=]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: tree format should be newick or nexus!\n"
	unless $tformat =~ /^(nexus|newick)$/i;
die " ERROR: provide a branch length cutoff if providing a tree!\n" 
	if $tree_in && ! $brlen_cut;
die " ERROR: provide a clades file!\n" unless $clades_in;


### MAIN
# loading clades #
my $clades_r = load_clades($clades_in);

# loading subclades #
my %subclades;
if($tree_in){
	my $treeo = tree_io($tree_in, $tformat);
	my $brlen_r = make_node_brlen_index($treeo);
		#print Dumper $brlen_r; exit;
	get_subclades($treeo, $brlen_r, $brlen_cut, \%subclades);
	}
else{		# each taxon is own subclade
	unique_subclades($clades_r, \%subclades);
	}

# making clade-subclade index #
my $cs_index_r = clade_subclade_index($clades_r, \%subclades);

# loading PA table
load_PA_table($clades_r, \%subclades, $cs_index_r, $subclade_cut, $absent_cut);
## grouping by subclade ##
## counting subclades ##


#-----Subroutines-----#
sub load_PA_table{
	my ($clades_r, $subclades_r, $cs_index_r, $subclade_cut, $absent_cut) = @_;

	# getting unique clades #
	my %uclades;
	map{$uclades{$_} = 1} values %$clades_r;
	print STDERR "Number of clades: ", scalar keys %uclades, "\n";
	
	# getting unique subclades #
	my %usubclades;
	map{$usubclades{$_} = 1} values %$subclades_r;
	print STDERR "Number of subclades: ", scalar keys %usubclades, "\n";
	
	# clade order for writing out table #
	my @clades = keys %$clades_r;
	
	# parsing PA table #
	my %header;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;

		# header parsing #
		if($.==1){
			for my $i (0..$#l){
				if($i > 2 && ! exists $clades_r->{$l[$i]}){		# should be taxa columns 
					print STDERR " WARNING: cannot find $l[$i] in clades list! $l[$i] not included!\n";
					next;
					}
				$header{$l[$i]} = $i;
				}
			# writing out header #
			print join("\t", @l[0..2], "fixed", @l[3..$#l]), "\n"; 	
				
			next;
			}
		
		# body #
		# determining if clades have ~all of gene #
		## grouping by subclade ##
		my %subclade_cnt;
		foreach my $taxon (keys %$subclades_r){
			next unless exists $header{$taxon};
			my $val = $l[$header{$taxon}];
			if($val =~ /NONE/){					# cluster not found in taxon
				$subclade_cnt{$subclades_r->{$taxon}}{'count'} += 0;
				}
			else{								# cluster found 
				$subclade_cnt{$subclades_r->{$taxon}}{'count'}++;
				}
			$subclade_cnt{$subclades_r->{$taxon}}{'total'}++;
			}
			
		## grouping by clade ##
		my %clade_cnt;
		foreach my $subclade (keys %subclade_cnt){
			# applying subclade cutoff #
			if($subclade_cnt{$subclade}{'count'} / $subclade_cnt{$subclade}{'total'}
				>= $subclade_cut){			# gene presence in >= fraction of subclde than cutoff
				$clade_cnt{$cs_index_r->{$subclade}}{'count'} += $subclade_cnt{$subclade}{'count'};
				}
			else{			
				$clade_cnt{$cs_index_r->{$subclade}}{'count'} += 0;
				}
			$clade_cnt{$cs_index_r->{$subclade}}{'total'}++;
			}
		
		# determining if cluster found in enough taxa to be called 'present in ~all members' #
		my %clade_PA;
		foreach my $clade (keys %clade_cnt){
			if($absent_cut >= 1){		# cutoff as number of abencenses allowed
				if($clade_cnt{$clade}{'total'} - $clade_cnt{$clade}{'count'} 
					<= $absent_cut){
					$clade_PA{$clade} = 1;
					}
				else{ $clade_PA{$clade} = 0; }
				}
			else{		# cutoff in fraction (number subclades present/total nunmber subclades)
				if( ($clade_cnt{$clade}{'total'} - $clade_cnt{$clade}{'count'}) / $clade_cnt{$clade}{'total'} 
					<= $absent_cut){		
					$clade_PA{$clade} = 1;
					}
				else{ $clade_PA{$clade} = 0; }		
				}
			}
		
		# determining if other clades have none of gene cluster #
		## counting by clade ##
		my %byclade_cnt;
		foreach my $taxon (keys %$clades_r){
			next unless exists $header{$taxon};
			my $val = $l[$header{$taxon}];
			if($val =~ /NONE/){					# cluster not found in taxon
				$byclade_cnt{$clades_r->{$taxon}}{'count'} += 0;
				}
			else{								# cluster found 
				$byclade_cnt{$clades_r->{$taxon}}{'count'}++;
				}
			}
		
		
		## getting number of clades w/ gene found ##
		my $n_clades_wGene;
		map{$n_clades_wGene++ if $byclade_cnt{$_}{'count'} > 0} keys %byclade_cnt;
		
				#print Dumper %clade_PA, "n_clades: $n_clades_wGene";
		
		# determining fixation #
		## only 1 clade w/ ~all taxa having cluster; other clades, must have none ##
		my %fixed;
		my $fixed = "NOT_FIXED";
		my $fixed_total = 0;
		foreach my $clade (keys %clade_PA){
			if($clade_PA{$clade} > 0 && $n_clades_wGene == 1){
				$fixed{$clade} = 1;
				$fixed = $clade;
				$fixed_total++;
				}
			else{ $fixed{$clade} = 0; }
			}
		die " ERROR: too many clades are fixed!\n"
			unless $fixed_total <2;

		# writing out line #
		print join("\t", @l[0..2], $fixed, @l[3..$#l]), "\n"; 
		}
	}

sub clade_subclade_index{
	my ($clades_r, $subclades_r) = @_;
	my %index;			# subclade_num => clade_ID
	#my %index;			# clade_id => subclade_num
	foreach my $taxon (keys %$clades_r){
		unless(exists $subclades_r->{$taxon}){
			print STDERR " WARNING: cannot find $taxon in subclade list! Skipping\n";
			next;
			}
			
		$index{$subclades_r->{$taxon}} = $clades_r->{$taxon};
		#$index{$clades_r->{$taxon}} = $subclades_r->{$taxon};
		}
		
	#	print Dumper %index; exit;
	return \%index;
	}

sub unique_subclades{
	my ($clades_r, $subclades_r) = @_;
	
	my $clade_cnt = 0;
	foreach my $taxon (keys %$clades_r){
		$clade_cnt++;
		$subclades_r->{$taxon} = $clade_cnt;
		}
	}

sub get_subclades{
# getting subclades defined by brlen cutoff #
	my ($treeo, $brlen_r, $brlen_cut, $subclades_r) = @_;
		
	# removing nodes from < height to > height #
	my $subclade_cnt = 0;
	#foreach my $node (sort{$brlen_r->{$a}<=>$brlen_r->{$b}} keys %$brlen_r){
	for my $node ($treeo->get_nodes){
		
		my $node_id;
		if( $node->is_Leaf){ $node_id = $node->id; }
		else{ $node_id = $node->description; }
		
		# removing taxa w/ height < brlen cutoff #
		unless($node->is_Leaf){
			next unless exists $brlen_r->{$node_id};
			next if $brlen_r->{$node_id} > $brlen_cut;		# node closer to tip then cutoff
			}
		
		#print Dumper $node_id;
		#print Dumper $node->is_Leaf;
		#print Dumper $nodes[0];
		
		# getting ancestor and seeing if it has brlen >= cutoff, if so calling this a subclade #
		my $anc = $node->ancestor;
		next unless exists $brlen_r->{$anc->description};
		next unless $brlen_r->{$anc->description} > $brlen_cut;	

		# subclade found; getting leaves #
		## all descendents in subclade ##
		$subclade_cnt++;
		if($node->is_Leaf){
			$subclades_r->{$node->id} = $subclade_cnt;
			}
		else{
			for my $child ( $node->get_all_Descendents ){
				next unless $child->is_Leaf;
				die " ERROR: ", $child->id, " already in subclade!\n"
					if exists $subclades_r->{$child->id};
				$subclades_r->{$child->id} = $subclade_cnt;
				}
			}
			
		}
	
	# printing clade distribution #
	print STDERR join("\t", qw/taxon subcladeID/), "\n";
	foreach my $taxon (sort {$subclades_r->{$a}<=>$subclades_r->{$b}} keys %$subclades_r){
		print STDERR join("\t", $taxon, $subclades_r->{$taxon}), "\n";
		}
	print STDERR "\nbug warning: root probably not included in a subclade!\n\n"; 
	
		#print Dumper %$subclades_r; exit;
	}

sub make_node_brlen_index{
# brlen index of internal nodes; brlen=max length from any leaf #
	my ($treeo) = @_;
	
	my %brlen;
	my $node_cnt = 0;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		# labeling internal nodes w/ unique ID #
		$node_cnt++;
		$node->description($node_cnt);
		next unless $node->ancestor;
		
		# getting distance from #
		my @brlens;
		for my $child ($node->get_all_Descendents){
			push @brlens, $treeo->distance(-nodes => [$node, $child]);
			}
		$brlen{$node->description} = max @brlens; #$treeo->distance(-nodes => [$node->ancestor, $node]);
		}

		#print Dumper %brlen; exit;
	return \%brlen;
	}
	
sub load_clades{
	my ($clades_in) = @_;
	
	open IN, $clades_in or die $!;
	my %clades;
	while(<IN>){
		chomp;
		my @l = split /\t/;
		die " ERROR: clades table does not have >=2 columns!\n"
			unless scalar @l >= 2;
			
		$clades{$l[0]} = $l[1]; 		# taxon=>clade
		}
	close IN;
	
	return \%clades;
	}

sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	$format =~ tr/[A-Z]/[a-z]/;
	my $input = Bio::TreeIO->new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		
	return $treeio;
	}

__END__

=pod

=head1 NAME

ITEP_PAfixed.pl -- determining gene clusters found in (~)all of 1 clade, but not in other clades

=head1 SYNOPSIS

ITEP_PAfixed.pl [flags]< ITEP_pres-abs.txt > ITEP_pres-abs_fixed.txt

=head2 required flags

=over

=item -clades  <char>

File designating clades (2-column; tab-delim; taxon\tclade).

=back 

=head2 options

=over

=item -tree  <char>

Tree file for determining subclades (newick or nexus).

=item -format  <char>

Tree file format (newick or nexus). [newick]

=item -brlen  <float>

Subclades defined as any clades where max branch length from any leaf to LCA is < '-brlen'.

=item -absent  <float>

Number of taxa in a clade that can be missing a gene in cluster (if < 1, fraction of cluster). [0]

=item -subclade  <float>

Fraction of taxa in subclade that must have cluster (eg. 0.51 for majority rules). [0]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_PAfixed.pl

=head1 DESCRIPTION

Add a column (4th column) which indicates gene cluster 'fixation' 
among clades. The cutoffs can help if the genomes are draft
and so some taxa may be artificially missing the gene cluster.

=head2 Fixation definition

The gene cluster must be found in enough taxa (after accounting
for absences in all taxa in clade or in subclades) to be considered
present throughout the clade.
Also, the cluster must be absent in all other clades.


=head1 EXAMPLES

=head2 No subclades defined:

ITEP_PAfixed.pl -c clades.txt < ITEP_pres-abs.txt > ITEP_pres-abs_fixed.txt

=head2 Defining subclades (collapsing tree w/ brlen cutoff of 0.001)

ITEP_PAfixed.pl -t tree.nwk -b 0.001 -c clades.txt < ITEP_pres-abs.txt > ITEP_pres-abs_fixed.txt

=head2 Majority rule for subclades

ITEP_PAfixed.pl -t tree.nwk -b 0.001 -s 0.51 -c clades.txt < ITEP_pres-abs.txt > ITEP_pres-abs_fixed.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

