#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use List::Util qw/min max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $prefix = "ranger-dtl_parse";
my $trans_cutoff = 0.025;
GetOptions(
	   "prefix=s" => \$prefix,
	   "transfer=f" => \$trans_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
# load table #
my ($trans_r, $dl_r, $Ntrees) = load_dtl_parsed();

# removing low-abundant transfers #
apply_trans_cutoff($trans_r, $trans_cutoff, $Ntrees);

# making ITOL tables #
make_hgt_table($trans_r, $prefix);
make_dl_table($dl_r, $prefix);

### Subroutines
sub make_dl_table{
# making duplication (& loss) for ITOL #
	my ($dl_r, $prefix) = @_;
	
	# summing #
	my @dup_vals;
	foreach my $node (keys %$dl_r){
		push @dup_vals, $dl_r->{$node}{"duplication"};
		}
	my $dup_max = max @dup_vals;


	# output FH #
	open OUT, ">$prefix\_Dup-Loss.txt", or die $!;

	# header #
	print OUT join("\t", qw/LABELS duplications/), "\n";
	print OUT join("\t", "COLORS", "#FF0000"), "\n";
	
	# writing each line #
	foreach my $node (keys %$dl_r){
		if($node =~ /\|/){ 		# internal node
			print OUT join("\t", $node, 
					join("", "R", 
						int($dl_r->{$node}{"duplication"} / $dup_max * 100)),
					$dl_r->{$node}{"duplication"}
					), "\n";
				}
		else{
			print OUT join("\t", $node, 
					$dl_r->{$node}{"duplication"}
					), "\n";
			}
		}
	close OUT
	}

sub make_hgt_table{
	my ($trans_r, $prefix) = @_;
	
	unless (%$trans_r){
		print STDERR " WARNING: no transfers above cutoff! Not transfer ITOL table written!\n";
		return 0;
		}
	
	# max/min values #
	my $min_trans = min values %$trans_r;
	my $max_trans = max values %$trans_r;
	
	# making color index #
	my $colors_r = heatmap_colors(10, $min_trans, $max_trans);

	# opening FH #
	open OUT, ">$prefix\_transfers.meta";

	# writing out HGT table #
	foreach my $trans (keys %$trans_r){
		foreach my $i (0..$#$colors_r){
			if($trans_r->{$trans} >= $$colors_r[$i][0]){
				print OUT join("\t", $trans, $$colors_r[$i][1]), "\n";
				}
			}
		}

	close OUT;
	}

sub apply_trans_cutoff{
# removing low abundant transfers #
## fraction of trees that predict transfer at a node ##
	my ($trans_r, $trans_cutoff, $Ntrees) = @_;
	print STDERR "Number of trees: $Ntrees\n";
	print STDERR "Total number of transfers (summed by node): ", scalar keys %$trans_r, "\n";
	
	#my $trans_sum = 0;
	#map{ $trans_sum += $_ } values %$trans_r;
	#my $trans_max = max values %$trans_r;
	
	my $del_cnt = 0;
	foreach my $trans (keys %$trans_r){
		if ($trans_r->{$trans} / $Ntrees < $trans_cutoff){
			delete $trans_r->{$trans};
			$del_cnt++;
			}
		}

	print STDERR "Number of tranfers (summed by node) below cutoff: $del_cnt\n";
	print STDERR "Number of transfers (summed by node) remaining: ", scalar keys %$trans_r, "\n";
		#print Dumper %$trans_r; exit;
	}

sub load_dtl_parsed{
# loading dtl_parsed table #
## (node -> node) => N-trees ##
## node => N-dupliations ##
	my %transfers;
	my %dl;
	my %trees; 
	while(<>){
		chomp;
		
		my @line = split /\t/;		# treeID, genetree_node, species_node, cat, recipient
		if($line[3] =~ /transfer/i){
			$transfers{join("\t", $line[2], $line[4])}++;
			}
		elsif($line[3] =~ /duplication/i){
			$dl{$line[2]}{"duplication"}++;
			}
			
		$trees{$line[0]} = 1;
		}
		print Dumper %transfers; exit;
	return \%transfers, \%dl, scalar keys %trees;
	}

sub heatmap_colors{
# color gradient of blue #
	my ($N, $min, $max) = @_;
	my @colors = qw/CCFFFF 99FFFF 33FFFF 0099CC 3366FF 0000FF 0000CC 000099 000066 000033/;
	map{$_ =~ s/^/#/} @colors;
	
	# span for colors #
	my @ret;
	for(my $i=0; $i<=9; $i+=10/$N){
		push @ret, int ($i + 0.5);		# rounding
		} 
	
	# color index #
	my @index;
	my $cnt = 0;
	for (my $i=$min; $i<=$max; $i+=($max-$min)/$N){
		last if $cnt == 9;
		push( @index, [ int ($i + 0.5), $colors[$cnt] ] );
		$cnt++;
		}

		#print Dumper @index; exit;	
	return \@index;
	}


__END__

=pod

=head1 NAME

ranger_dtl_parse_ITOL.pl -- convert ranger-dtl parse output into ITOL metadata

=head1 SYNOPSIS

ranger_dtl_parse_ITOL.pl < ranger_dtl_parsed.txt

=head2 options

=over

=item -p 	Output file prefix. ["ranger-dtl_parse"]

=item -t 	Cutoff for fraction of trees predicting transfer at a node. [0.025]

=item -h	This help message

=back

=head2 For more information:

perldoc ranger_dtl_parse.pl

=head1 DESCRIPTION

Make an ITOL metadata table of duplications and also
an ITOL table of horzontal gene transfers (HGTs).
The HGT table will be colored by blue gradient from 
min number of trees supporting the HGT to the max
(range is defined after applying the cutoff, so
it is just the range of the HGTs plotted).

For the duplication-loss table, the radius ('R#')
column is defined relative to the node with the
max number of trees supporting the duplication.

=head1 EXAMPLES

=head2 Usage: 

ranger_dtl_parse_ITOL.pl < ranger-dtl_parsed.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

