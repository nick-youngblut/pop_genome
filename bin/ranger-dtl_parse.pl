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

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
parse_ranger_dtl();

### Subroutines
sub parse_ranger_dtl{
# parsing output from ranger-dtl #
	my %dtl;
	my %lca;
	my $tree_id;
	while(<>){
		chomp;
		next if /^\s*$/;

		if(/ Reconciliation for Gene Tree/){
			($tree_id = $_) =~ s/.+Reconciliation for Gene Tree (\d+).+/$1/;
			read_species_tree(\%lca) unless %lca;
			}
		elsif(/^Reconciliation:/){
			parse_dtl_block(\%dtl, \%lca, $tree_id);
			}
		}
	}
	
sub parse_dtl_block{
# parsing DTL portion of output #
	my ($dtl_r, $lca_r, $tree_id) = @_;
	
	while(<>){
		chomp;
		last if /^\s*$/;
		
		if(/Leaf Node/){		# if leaf in gene tree
			(my $line = $_) =~ s/: .+//;
			#$dtl_r->{$tree_id}{"Leaf"} = [($line)];
			print join("\t", $tree_id, $line, "", "Leaf", ""), "\n";
			}
		else{					# if internal gene tree node
			my @line = split/= +LCA\[|\]: +|, +| +--> +/;
			
			# mapping to LCA or leaf? #
			my $mapping;
			if( exists $lca_r->{$line[5]} ){ 
				$mapping = join("|", @{$lca_r->{$line[5]}});
				}
			else{ $mapping = $line[5]; }
			
			# recipient LCA or leaf? #
			my $recip;
			if( $line[7] && exists $lca_r->{$line[7]}){
				$recip = join("|", @{$lca_r->{$line[7]}} );
				}
			else{ $recip = $line[7]; }

			# writing out LCA lines #
			if(/Speciation/){
				print join("\t", $tree_id, join("|", @line[1..2]),
							$mapping, $line[3], ""), "\n";
				}
			elsif(/Duplication/){
				print join("\t", $tree_id, join("|", @line[1..2]),
							$mapping, $line[3], ""), "\n";
				}
			elsif(/Transfer/){
				print join("\t", $tree_id, join("|", @line[1..2]),
							$mapping, $line[3], $recip), "\n";			
				}
			}
		}
	}

sub read_species_tree{
# reading species tree to making node_ID -> LCA hash #
	my ($lca_r) = @_;
	
	# loading species tree as string #
	my $tree_str;
	while(<>){
		if(/^Species Tree:/){
			while(<>){
				$tree_str .= $_;
				last if /;\s*$/;
				}
			}
		last if $tree_str;
		}
	
	# getting LCA for each node using bioperl #
	open (my $str_fh, "<", \$tree_str) or die $!;
	my $treeio = Bio::TreeIO->new( -fh => $str_fh, -format => "newick");

	my $node_cnt;
	while (my $tree = $treeio->next_tree){
		
		for my $node ($tree->get_nodes){
			next if $node->is_Leaf;
			my $node_id = $node->id;
			$node_cnt++;

			# getting all leaves of node #
			my @children;
			for my $child ($node->get_all_Descendents){
				push @children, $child if $child->is_Leaf;
				}
			
			# finding a pair of leaves where LCA is node #
			for (my $i=0; $i<=$#children; $i++){
				last if exists $lca_r->{$node_id};
				for (my $ii=$#children; $i>=0; $i--){			# working backward
					last if exists $lca_r->{$node_id};
					
					my $lca_id = $tree->get_lca( @children[($i, $ii)] )->id;
					if($lca_id =~ /[A-Za-z]/ && $lca_id eq $node_id){
						$lca_r->{$node_id} = [($children[$i]->id, $children[$ii]->id)];
						}
					elsif( $lca_id == $node_id){
						$lca_r->{$node_id} = [($children[$i]->id, $children[$ii]->id)];
						}
					else{ die $!; }
					}
				}
			
			$lca_r->{$node_id} = "NA" unless exists $lca_r->{$node_id};
			}
		}

	# sanity check #
	die " ERROR: LCAs not found for all internal nodes in species tree!\nWere bootstrap values left in the trees???\n"
		unless $node_cnt == scalar keys %$lca_r;
	
		#print Dumper %$lca_r; exit;
	return $lca_r;
	}


__END__

=pod

=head1 NAME

ranger_dtl_parse.pl -- parse ranger-dtl output into a *txt table

=head1 SYNOPSIS

ranger_dtl_parse.pl < ranger.dtl.out > ranger.dtl.out.txt

=head2 options

=over

=item -h	This help message

=back

=head2 For more information:

perldoc ranger_dtl_parse.pl

=head1 DESCRIPTION

The output is a tab-delimited table. Blank value = not applicable.

=head2 Columns of output

=over

=item * 	Tree ID (order in newick input to ranger-dtl)

=item * 	Gene tree node name

=item * 	Species tree node name

=item * 	Category (Duplication, Transfer, Speciation, or Leaf)

=item * 	Recipient (species tree node name)

=back

=head1 EXAMPLES

=head2 Usage: 

ranger_dtl_parse.pl < ranger-dtl.out > ranger-dtl.out.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

