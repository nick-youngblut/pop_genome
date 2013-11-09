#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $tree_in, $norm_bool);
my $format = "newick";
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "normalize" => \$norm_bool, 		# normalize counts by cluster
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if($tree_in){
	die " ERROR: $tree_in file not found!\n" unless -e $tree_in;
	$format = format_check($format);
	}


### MAIN
# loading gene info table #
my ($geneInfo_r, $clusters_r) = load_gene_info();

# making basic metadata table #
my $cols_r = get_hex_colors(scalar @$clusters_r);
normalize_counts($geneInfo_r, $clusters_r) if $norm_bool;
write_leaf_meta($geneInfo_r, $clusters_r, $cols_r);

# adding internal node metadata (sums of external nodes) #
if($tree_in){
	my $treeo = tree_io($tree_in) if $tree_in; # loading tree #
	# unique node ids #
	unique_node_IDs($treeo);
	
	# getting counts per node & writing #
	get_node_counts($treeo, $geneInfo_r, $clusters_r);
	}

# normalizing each org by total COG #
#normalize_COG(\%COG);

# writing out metadata table #
#write_metadata(\%COG, $labels_r);


### Subroutines
sub get_node_counts{
# getting counts per node & writing #
	my ($treeo, $geneInfo_r, $clusters_r) = @_;
	
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		
		my @children; 
		for my $child ($node->get_all_Descendents){
			next unless $child->is_Leaf;
			push @children, $child;
			}
			
		# finding leaf pair w/ correct LCA #
		my $node_id = $node->id;
		my @leaves;
		my $last_bool = 0;
		for my $i (0..$#children){
			last if $last_bool;
			for my $ii (0..$#children){
				next if $i >= $ii;		# all pairwise (lower triange)
				my $lca = $treeo->get_lca(-nodes => [$children[$i], $children[$ii]]);
				if($lca->id eq $node_id){
					push @leaves, $children[$i]->id, $children[$ii]->id;
					$last_bool = 1;
					last;
					}
				}
			}
		
		# getting totals for all leaves #
		my %counts;
		foreach my $child (@children){
			my $child_id = $child->id;
			foreach my $cluster (@$clusters_r){
				if(exists $geneInfo_r->{$child_id}{$cluster}){
					#print Dumper $geneInfo_r->{$child_id}{$cluster}; exit;
					$counts{$cluster} += $geneInfo_r->{$child_id}{$cluster};
					$counts{'total'} += $geneInfo_r->{$child_id}{$cluster};
					}
				else{
					$counts{$cluster} += 0;
					$counts{'total'} += 0;
					}
				}
			}
			
		# writing out line #
		#$counts{'total'} = 1 unless $counts{'total'};
		print join("\t", join("|", @leaves), "R$counts{'total'}", 
					@counts{@$clusters_r}), "\n";
		}

	}
	
sub unique_node_IDs{
	my ($treeo) = @_;
	
	my $node_cnt = 0;
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		$node_cnt++;
		$node->id($node_cnt);
		}
	}

sub write_leaf_meta{
# writing counts for each cluster #
	my ($geneInfo_r, $clusters_r, $cols_r) = @_;
	
	# header #
	print join("\t", "LABELS", @$clusters_r), "\n";
	print join("\t", "COLORS", @$cols_r), "\n";
	
	# body #
	foreach my $taxon (keys %$geneInfo_r){
		my @l;
		foreach my $cluster (@$clusters_r){
			if(exists $geneInfo_r->{$taxon}{$cluster}){
				push @l, $geneInfo_r->{$taxon}{$cluster};
				}
			else{ push @l, 0; }
			}
		print join("\t", $taxon, @l), "\n";
		}
	}

sub normalize_counts{
# normalizing by cluster count total # 
	my ($geneInfo_r, $clusters_r) = @_;

	foreach my $cluster (@$clusters_r){	
		# getting total #
		my $total = 0;
		foreach my $taxon (keys %$geneInfo_r){
			$total += $geneInfo_r->{$taxon}{$cluster} 
				if exists $geneInfo_r->{$taxon}{$cluster};
			}
		# normalizing by cluster #
		foreach my $taxon (keys %$geneInfo_r){
			$geneInfo_r->{$taxon}{$cluster} = $geneInfo_r->{$taxon}{$cluster} / $total
				if exists $geneInfo_r->{$taxon}{$cluster};
			}		
		}

	#print Dumper %$geneInfo_r; exit;
	}

sub load_gene_info{
# loading gene info table #
	my %geneInfo;
	my %gene_clusters;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die " ERROR: line $. doesn't have 14 columns!\n"
			unless scalar @l >= 14;
		
		$geneInfo{$l[1]}{$l[13]}++;
		$gene_clusters{$l[13]} = $. unless exists $gene_clusters{$l[13]};
		}
	
		#print Dumper %geneInfo; exit;
		#print Dumper %gene_clusters; exit;
	return \%geneInfo, [sort {$gene_clusters{$a}<=>$gene_clusters{$b}} keys %gene_clusters];
	}

sub get_hex_colors{
# getting hexidecimal colors for metadata #
	my $n_col = shift;
	$n_col--;
	my @hex = qw/FF0000 FF6600 33FF00 0000FF FF00FF FF0099 33CCFF 990000 CC0099 000066 006600 CC6600/;
	map{$_ =~ s/^/#/} @hex;
	
	if($n_col > scalar @hex){
		for my $i (0..int( $n_col/ scalar @hex)){
			push @hex, @hex;
			}
		}
	
	return [@hex[0..$n_col]];
	}


sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO->new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;		
	return $treeio;
	}

sub format_check{
	my $format = shift;
	if ($format =~ /nwk|newick/i){ $format = "newick" if $format =~ /nwk|newick/i; }
	elsif ($format =~ /nex|nexus|tre/i){ $format = "nexus"; }
	else{ die " ERROR: format not recognized [newick or nexus accepted]\n"; }
	return $format;
	}

__END__

=pod

=head1 NAME

ITEP_geneInfo2ITOL.pl -- ITEP gene info table to an ITEP metadata table

=head1 SYNOPSIS

ITEP_geneInfo2ITOL.pl [options] < ITEPgeneInfo.txt > ITOL.meta

=head2 Options

=over

=item -tree  <char>

Tree file (newick or nexus).

=item -format  <char> 	

Tree format (newick or nexus). [newick] 

=item -normalize  <bool>

Normalize by cluster? [FALSE]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfo2ITOL.pl

=head1 DESCRIPTION

Make an ITOL metadata table of the number of genes in a cluster
for each taxon. If a tree is provided, leaf counts are totaled
for each node and added to the metadata table for displaying pie charts.

=head1 EXAMPLES

=head2 Basic usage:

ITEP_geneInfo2ITOL.pl < ITEPgeneInfo.txt > ITOL.meta

=head2 Normalizing by cluster:

ITEP_geneInfo2ITOL.pl -n < ITEPgeneInfo.txt > ITOL.meta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

