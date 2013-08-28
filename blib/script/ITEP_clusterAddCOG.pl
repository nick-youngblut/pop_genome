#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $peg_col = 1;
GetOptions(
	   "fig=i" => \$peg_col,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$peg_col--;

### MAIN
while(<>){
	chomp;
	next if /^\s*$/;
	my @line = split /\t/;
	
	die " ERROR: PEG column not found!\n"
		unless $line[$peg_col];
	my $pegID = $line[$peg_col];
	
	my $cog_cat = get_external_hit($pegID);
	print join("\t", @line, $cog_cat), "\n";
	}

#--- old ---#
#my $cluster_r = load_cluster_table($ARGV[0], $peg_col);
#my $hit_index_r = get_external_hit_index($ARGV[0], $peg_col, $cluster_r);
#my $cog_r = get_cog_info($ARGV[0], $peg_col, $cluster_r);

#join_tables($cluster_r, $hit_index_r, $cog_r);

#-----Subroutines-----#
sub get_external_hit{
	my ($pegID) = @_;
	
	my $cmd = "echo \"$pegID\" | db_getExternalClusterGroups.py -d COG | \
	db_getExternalClustersById.py -c 13 ";
	
	my $out = `$cmd`;
	my @out = split /\n|\r/, $out;
	
	# parsing hit #
	my %cog_cat;
	foreach (@out){
		my @line = split /\t/;
		(my $COG_cat = $line[3]) =~ s/.*\[|\].*//g;
		$cog_cat{$COG_cat}++;
		}
	
	if(%cog_cat){ return join(";", keys %cog_cat); }
	else{ return ""; }
	}

sub join_tables{
# 3 table join #
	my ($cluster_r, $hit_index_r, $cog_r) = @_;
	
	# making fig=>cdd_id=>cog
	my %index;
	foreach my $fig (keys %$hit_index_r){
		foreach my $cdd_id ( keys %{$hit_index_r->{$fig}} ){
			$index{$fig}{$cdd_id} = $cog_r->{$cdd_id}
				if exists $cog_r->{$cdd_id}; 
			}
		}
	
	# adding to original table #
	foreach my $fig (keys %$cluster_r){
		next unless exists $index{$fig};
		
		foreach my $row (@{$cluster_r->{$fig}}){
			foreach my $cdd_id (keys %{$index{$fig}}){
				foreach my $cog_id (keys %{$index{$fig}{$cdd_id}}){
					print join("\t", @$row, @{$index{$fig}{$cdd_id}{$cog_id}}), "\n";
					}
				}
			}
		}
	}

sub get_cog_info{
	my ($infile, $peg_col, $cluster_r) = @_;
	
	$peg_col++;
	
	open PIPE, "cat $infile | cut -f $peg_col | db_getExternalClusterGroups.py | db_getExternalClustersById.py -c 13 |" or die $!; #| 
	
	my %cog;
	while(<PIPE>){
		chomp;
    	my @line = split /\t/;
 			#print Dumper @line; exit;   	
    	next unless $line[1] =~ /COG\d+/;
    	my @cat = split /\[/, $line[3];

    	map{s/\]//g} @cat;    	
    	push(@line, @cat);
    	
    	#push(@{$cog{$line[0]}}, [@line[1..2,5..6]]);			# fig=>cdd_id
		$cog{$line[0]}{$line[1]} = \@cat;
		}
	close PIPE;	

		#print Dumper %cog; exit;	
	return \%cog;
	}

sub get_external_hit_index{
	my ($infile, $peg_col, $cluster_r) = @_;
	
	$peg_col++;
	
	open PIPE, "cat $infile | cut -f $peg_col | db_getExternalClusterGroups.py |" or die $!; #| db_getExternalClustersById.py -c 13
	
	my %hit_index;
	while(<PIPE>){
		chomp;
    	my @line = split /\t/;
 			#print Dumper @line; exit;   	
    	$hit_index{$line[0]}{$line[1]} = 1;			# fig=>cdd_id
		}
	close PIPE;	

#		print Dumper %hit_index; exit;	
	return \%hit_index;
	}

sub load_cluster_table{
# load gene cluster table #
	my ($infile, $peg_col) = @_;
	
	open IN, $infile or die $!;
	
	my %cluster;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @line = split /\t/;
		
		push( @{$cluster{$line[$peg_col]}}, \@line);
		}
	close IN;
		#print Dumper %cluster; exit;
	return \%cluster;
	}



__END__

=pod

=head1 NAME

ITEP_clusterAddCOG.pl -- adding COG info to table containing PEG IDs (e.g. geneInfoTable)

=head1 SYNOPSIS

ITEP_clusterAddCOG.pl [options] < pegIDs.txt > pegIDs_COG.txt

=head2 options

=over

=item -p 	PEG (geneID) column (indexed by 1). [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_clusterAddCOG.pl

=head1 DESCRIPTION

Use ITEP to appending COG functional categories onto the end of the
provided list/table.

Uses db_getExternalClusterGroups.py & db_getExternalClustersById.py

=head1 EXAMPLES

=head2 Basic usage:

ITEP_clusterAddCOG.pl < geneInfoTable.txt > geneInfoTable_COG.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

