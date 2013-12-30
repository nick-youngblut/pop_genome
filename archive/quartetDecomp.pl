#!/usr/bin/env perl

=pod

=head1 NAME

quartetDecomp.pl -- decompose all input trees into quartets

=head1 SYNOPSIS

quartetDecomp.pl [flags] > quartet_summary.txt

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

perldoc quartetDecomp.pl

=head1 DESCRIPTION

summarize output from Quartet Decomposition Server analysis

=head1 EXAMPLES

=head2 Basic usage:

quartetDecomp.pl [flags] > quartet_summary.txt

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
use Bio::TreeIO;
use List::Util qw/min/;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $trees_in, $boot_id_b);
my $boot_cut = 70;
GetOptions(
	"trees=s" => \$trees_in,
	"id" => \$boot_id_b, 
	"bootstrap=i" => \$boot_cut,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
die "ERROR: provide a file listing tree file names!\n"
	unless $trees_in;
die "ERROR: cannot find tree file list!\n"
	unless -e $trees_in;

#--- MAIN ---#
my $files_r = load_tree_list($trees_in);

foreach my $file (@$files_r){
	my $treeo = tree_io($file, 'newick', $boot_id_b);
	
	# adding '0' to node IDs if not found #
	add_boot($treeo);
	
	# all possible quartets #
	my ($sets_r, $index_r) = count_possible_quartets($treeo, $file);

	foreach my $set (@$sets_r){
		# rerooting tree for each taxon #
		my @dists;
		for my $i (0..$#$set){
			$treeo->reroot( $index_r->{$$set[$i]} ) 
				unless $treeo->get_root_node->id eq $index_r->{$$set[$i]}->id;
		
			# getting topographic distance to root for other 3 in set #
			for my $ii (0..3){
				next if $i == $ii;
				push @dists, [$i, $ii, 
					get_topo_dist_root($treeo, $index_r->{$$set[$ii]}, $index_r->{$$set[$i]}) ];
				}
			}
			
		print Dumper @dists; exit;
		# find min dist between any 2 taxa #
			# if >2 have min dist; multiple quartets needed 
			# other 2 taxa are other node
		# foreach pair in quartet: get lca bootstrap #
		
		# write out quartet w/ boostrap support values #
		}
	
	}

#--- Subroutines ---#
sub get_min_dist_node{
	my ($dists_r) = @_;
	
	my %cnt;
	map{ $cnt{$_}++ } @$dists_r;
	my $min_dist = min @$dists_r;
	
	#print Dumper @$dists_r; 
	if ($cnt{$min_dist} > 1){
		warn "warn: >=2 nodes have the same minimum topo distance from root! Skipping\n"
			if $verbose_b;
		return 0;
		}

	for my $i (0..$#$dists_r){
		return $i if $i == $min_dist;
		}
	}

sub get_topo_dist_root{
# number of interveining nodes between query taxon & root #
	my ($treeo, $query, $root) = @_;
	
	# labeling root ancestor ID #
	my $root_anc_id = $root->ancestor->id;
	$root->ancestor->id('root_anc');
	
	my $dist = 0;
	my $anc = $query;
	while(1){
		die "ERROR: no ancestor for ", $anc->id, "!\n"
			unless defined $anc->ancestor;
		last if $anc->ancestor->id eq $root->ancestor->id;
		$dist++;
		$anc = $anc->ancestor;
		$anc->id('None') unless defined $anc->id;
		}
		
	$root->ancestor->id($root_anc_id);
	
		#print Dumper $dist;
	return $dist;
	}

sub add_boot{
	my ($treeo) = @_;
	
	for my $node ($treeo->get_nodes){
		next if $node->is_Leaf;
		$node->id(0) unless defined $node->id;
		}
	}

sub count_possible_quartets{
# all sets of 4 taxa #
	my ($treeo, $file) = @_;
	
	my @taxa = $treeo->get_leaf_nodes;
	my %index;
	for my $i (0..$#taxa){
		$index{$i} = $taxa[$i];
		}
		
	my %u_tax;
	foreach my $i (0..$#taxa){
		foreach my $ii (0..$#taxa){
			foreach my $iii (0..$#taxa){
				foreach my $iv (0..$#taxa){
					
					my @tax = sort{$a<=>$b} ($i, $ii, $iii, $iv) ;
					$u_tax{join(" ", @tax)} = 1;
					}
				}
			}
		}
	
	print STDERR "Number of possible quartets for $file: ",
			scalar keys %u_tax, "\n" if $verbose_b;
	
	# array of all taxa sets #
	my @sets;
	map{ push @sets, [split / /, $_] } keys %u_tax;
	
		#print Dumper @sets; exit;
	return \@sets, \%index;
	}

sub tree_io{
	# loading tree object: just 1st tree #
	my ($tree_in, $format, $boot_id_b) = @_;
	my $input;
	if($boot_id_b){
		$input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
		}
	else{
		$input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format,
								-internal_node_id => 'bootstrap');
		}			
					
	my $treeo = $input->next_tree;	
	return $treeo;
	}

sub load_tree_list{
	my ($trees_in) = @_;

	open IN, $trees_in or die $!;
	my @files;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		die "ERROR: cannot find '$_'!\n" unless -e $_;
		push @files, $_;
		}
	return \@files;
	}
