#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults


### MAIN
# loading gene info table #
my ($geneInfo_r, $index_r) = load_gene_info();

# making basic metadata table #
#my $cols_r = get_hex_colors(scalar keys %{$index_r->{'cluster_order'}});

# making table 
writeTable($geneInfo_r, $index_r);



#--- Subroutines ---#

=head2 writeTable

Writing out tab-delimited table.
Columns: taxon, copy, cluster

=cut

sub writeTable{
  my ($geneInfo_r, $index_r) = @_;

  # header 
  my @clusters_ordered = sort{$index_r->{cluster_order}{$a} <=> 
				$index_r->{cluster_order}{$b}}
    keys %{$index_r->{cluster_order}};

  print join("\t", 'Taxon', 'Copy_number', @clusters_ordered), "\n";

  # body
  foreach my $taxon (keys %$geneInfo_r){
    foreach my $copy (keys %{$geneInfo_r->{$taxon}}){  # don't actually need max copy
      # NA if not found
      map{ $geneInfo_r->{$taxon}{$copy}{$_}{'gene_len'} = 'NA' unless exists
	     $geneInfo_r->{$taxon}{$copy}{$_}{'gene_len'} } @clusters_ordered;
      
      # array of gene lengths
      my @lengths =  map{ $geneInfo_r->{$taxon}{$copy}{$_}{'gene_len'} } @clusters_ordered;


      # writing
      print join("\t", $taxon, $copy, @lengths), "\n";
    }
  }
}


=head2 load_gene_info

=head3 IN

ITEP gene info via STDIN

=head3 OUT

geneInfo : {taxon : cluster : copy : ('gene_len'|'annotation') : value }
index :  {cluster : (order|copy) : value}

=cut

sub load_gene_info{
# loading gene info table #
  my %geneInfo;
  my %index;
  while(<>){
    chomp;
    next if /^\s*$/;
    my @l = split /\t/;
    die " ERROR: line $. doesn't have 14 columns!\n"
      unless scalar @l >= 14;

    # columns needed: 
    ## 2 = taxon_name
    ## 6 = start
    ## 7 = end
    ## 10 = annotation
    ## 14 = cluster
    
    my $taxonID = $l[1];
    my $gene_start = $l[5];
    my $gene_end = $l[6];
    my $annotation = $l[9];
    my $clusterID = $l[13];
    my $gene_len = abs($gene_end - $gene_start);

    # copy
    my $copy = 1;
    while(1){
      if(exists $geneInfo{$taxonID}{$copy}{$clusterID}){
	$copy++;
      }
      else{ last; }
    }
    
    # loading hash
    $geneInfo{$taxonID}{$copy}{$clusterID}{'gene_len'} = $gene_len;

    # all unique clusters
    $index{'unique_clusters'}{$clusterID} = 1;
    # copies per cluster
    $index{'copy'}{$clusterID}{$copy} = 1;
  }  

  # setting index: 
  ## ordering of clusterIDs (smallest to largest)
  my $cnt = 0;
  foreach my $clusterID (sort{$a<=>$b} keys %{$index{'unique_clusters'}}){
    $index{'cluster_order'}{$clusterID} = $cnt;    
    $cnt++;
  }
  ## max copy per cluster
  foreach my $clusterID (keys %{$index{'unique_clusters'}}){
    $index{'max_copy'}{$clusterID} = max keys %{$index{'copy'}{$clusterID}};
  }
  
#  print Dumper %index; exit;
#  print Dumper %geneInfo; exit;
  return \%geneInfo, \%index;
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



__END__

=pod

=head1 NAME

ITEP_geneInfo2GeneLen.pl -- ITEP gene info table to a table of gene lengths.

=head1 SYNOPSIS

ITEP_geneInfo2GeneLen.pl [options] < ITEPgeneInfo.txt > geneLengths.txt

=head2 Options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_geneInfo2GeneLen.pl

=head1 DESCRIPTION

Convert an ITEP gene info table into a
tab-delimited table of gene lengths.

Multi-copy genes will each be on a seperate line.
The grouping of multi-copy genes on each line is random.

=head2 OUTPUT

columsn: taxon_name, gene_copy_number, clusterID1, clusterID2, ...

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

