#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

ITEP_traceGeneClust.pl -- tracing gene-homolog relationships across clust runs

=head1 VERSION

This is version 0.0.1

=head1 USAGE

ITEP_traceGeneClust.pl [options]

=head1 REQUIRED ARGUMENTS

=item -g[eneList] <gl>

List of FIG-PEG IDs. '-' if from STDIN

=for Euclid:
gl.type: input


=item -c[[luster][Runs]] <cr>

comma-separated list of cluster runs to trace relationship.

=for Euclid:
cr.type: string


=head1 OPTIONS

=over

=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --quiet

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

Trace homology for a gene across ITEP cluster runs.
Homology tracing: geneID --> clusterID_run1 --> geneIDs_in_Cluster 
--(foreach clusterID)--> clusterID_run2 --> ...

This is need when a gene is not found in all cluster runs
in the homology trace. Thus, this will associate a geneID
with all clusters (in user-listed cluster runs) that have
homology via inclusion of any of the same genes.

The output is a tab-delimited table written to STDOUT.

=head2 Warnings:

ITEP must be in PATH.

=head1 AUTHOR

Nick Youngblut (ndy2@cornell.edu)

=head1 BUGS

There are undoubtedly serious bugs lurking somewhere in this code.
Bug reports and other feedback are most welcome.

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

use Data::Dumper;
use Getopt::Euclid;
use IPC::Cmd qw/can_run run/;


#--- I/O error ---#
my @clustRuns = split /,/, $ARGV{-clusterRuns};

map{ die "ERROR: cannot find '$_'\n" unless can_run($_) }
  qw/db_getClustersContainingGenes.py 
     db_getClusterGeneInformation.py/;

#--- MAIN ---#
# workflow:
#  foreach geneID:
#   recursively get clusters for cluster

# getting gene list
my $geneList_r = readGeneList($ARGV{-geneList});

# calling ITEP to recursively get clusters 
my %traces;
foreach my $geneID (@$geneList_r){
  $traces{$geneID} = getCluster($geneID, undef, \@clustRuns);
}

# writing trace table
print join("\t", 'geneID', @clustRuns), "\n";
writeTrace(\%traces);


#--- Subroutines ---#
sub writeTrace{
  # recursively writing trace for each geneID
  my $traces_r = shift;

  foreach my $geneID (keys %$traces_r){
    getTrace($traces_r->{$geneID}, $geneID);
  }
}


sub getTrace{
  my $clusthash_r = shift;
  my $geneID = shift;
  my $path_r = shift;
  

  if( ref($clusthash_r) eq 'HASH' ){    
    foreach my $clustID (keys %$clusthash_r){
      push @$path_r, $clustID;
      getTrace($clusthash_r->{$clustID}, $geneID, $path_r);
    }
  }
  else{ # end of  path; writing
    print join("\t", $geneID, @$path_r), "\n";
  }
}


sub getTraceLevels{
  my $clusthash_r = shift;
  my $level = shift;

  $level = 0 unless defined $level;
  $level++;

#  my $paths_r = shift;

  if( ref($clusthash_r) eq 'HASH' ){    
    foreach my $clustID (keys %$clusthash_r){
      print "$level\t$clustID\n";
      getTrace($clusthash_r->{$clustID}, $level);
    }
  }
  else{
    print "$level\tEND\n";
#    return $paths_r;
  }
}


sub getCluster{
  # recursively tracing clusters 
  my $geneID = shift;
  my $clustID = shift;
  my $clustRuns_r = shift;
  my $clustRunIdx = shift;

  # input check
  $clustRunIdx = 0 unless defined $clustRunIdx;
  my $nextClustRun = defined $clustRuns_r->[$clustRunIdx] ?
    $clustRuns_r->[$clustRunIdx] :
      return 0;

  # cmd
  my $cmd;
  if(defined $geneID){
    $cmd = "echo '$geneID' | db_getClustersContainingGenes.py | grep $nextClustRun |";
  }
  elsif(defined $clustID){
    my $lastClustRun = $clustRuns_r->[$clustRunIdx -1];    
    $cmd = "printf '$lastClustRun\\t$clustID' | db_getClusterGeneInformation.py | \
            cut -f 1 | db_getClustersContainingGenes.py | grep $nextClustRun |";
  }
  else{
    die "ERROR: geneID or clustID must be defined!\n";
  }

  # parsing output
  open PIPE, $cmd or die $!;
  my %uclusts;
  while(<PIPE>){
    chomp;
    next if /^\s*$/;
    
    my @l = split /\t/;
    die "ERROR: line $. is not 3 columns!\n"
      unless scalar @l == 3;
    
    $uclusts{$l[1]} = 1;
  }
  close PIPE;

  # saving and recursively calling
  $clustRunIdx++;
  my %clust2clust;
  foreach my $clustID (keys %uclusts){
    $clust2clust{$clustID} = getCluster(undef, $clustID, 
					$clustRuns_r, $clustRunIdx);
  }

  return \%clust2clust;
}


sub readGeneList{
  my $infile = shift;
  
  my $ifh;
  $infile eq '-' ? $ifh = \*STDIN :
    open $ifh, $infile or die $!;

  my @geneList;
  while(<$ifh>){
    chomp;
    next if /^\s*$/;
    
    my @l = split /\t/;
    push @geneList, $l[0];
  }
  close $ifh or die $!;

  return \@geneList;
}
