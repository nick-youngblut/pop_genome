#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Parallel::ForkManager;
use Bio::TreeIO;
use IPC::Cmd qw/can_run/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $rename_b, $clusters_in);
my $threads = 2;
my $fork = 0;
my $prefix = "clusters";
my $raxml_prog = "raxmlHPC-PTHREADS-SSE3";
my $mafft_prog = "mafft-linsi";
my $pal2nal = "pal2nal.pl";
GetOptions(
	   "rename" => \$rename_b, 			# call replaceOrgWithAbbrev? [TRUE]
	   "threads=i" => \$threads, 		# number of threads to cal
	   "fork=i" => \$fork, 				# number of parallel files to process
	   "prefix=s" => \$prefix, 			# output dir prefix
	   "clusters=s" => \$clusters_in, 	# cluster list
	   "raxml=s" => \$raxml_prog, 		# which raxml to use?
	   "mafft=s" => \$mafft_prog, 		# which mafft to use?
	   "pal2nal=s" => \$pal2nal, 		# name pal2nal executable
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a list of clusters!\n"
	unless $clusters_in;
die " ERROR: cannot find $clusters_in!\n"
	unless -e $clusters_in;
my $curdir = File::Spec->rel2abs(File::Spec->curdir());

# checking executables #
can_run($pal2nal) or die "ERROR: Cannot call '$pal2nal'. Is it in your PATH?\n";
can_run($raxml_prog) or die "ERROR: Cannot call '$raxml_prog'. Is it in your PATH?\n";
can_run($mafft_prog) or die "ERROR: Cannot call '$mafft_prog'. Is it in your PATH?\n";

### Main
# calling ITEP scripts to get fastas #
get_fastas($clusters_in, $curdir, $prefix);

# making directories # 
my ($align_dir, $pal2nal_dir, $ML_dir, $rn_dir) = make_dirs($curdir, $prefix);

# getting list of fasta files #
my $fasta_r = get_fasta_files($curdir, $prefix);

# forking #
my $pm = new Parallel::ForkManager($fork);
foreach my $cluster (@$fasta_r){
  print STDERR "Starting alignment & phylogeny on: '$cluster'\n" unless $verbose;
  $pm->start and next;	# starting fork
  # purging fasta names of RAxML-unfriendly characters #
  purge_names($cluster, $curdir);
  
  # changing any internal stops in AA sequences '*' to 'X'
  edit_stops($cluster, $curdir);
  
  # getting number of sequences in file #
  my $seq_cnt = get_seq_cnt($cluster, $curdir);
  
  # checking for >= 2 taxa in cluster #
  die "WARNING: <2 sequences in '$cluster'. Exiting prior to calling mafft!\n"
    unless $seq_cnt >= 2;
  
  # AA alignment #
  call_mafft($cluster, $mafft_prog, $curdir, $prefix, $align_dir, $threads);
  # pal2nal #
  call_pal2nal($cluster, $curdir, $prefix, $align_dir, $pal2nal_dir, $rename_b, $pal2nal);	
  
  # checking for >= 4 taxa in cluster (required for RAxML) #
  die "WARNING: <4 sequences in '$cluster'. Exiting prior to calling RAxML!\n"
    unless $seq_cnt >= 4;
  
  # phylip & raxml #
  phy2raxml($cluster, $curdir, $prefix, $pal2nal_dir, $ML_dir, $threads, $raxml_prog);
  # removing PEGs (for comparing to species tree) #
  rm_PEGs($cluster, $curdir, $prefix, $ML_dir, $rn_dir);
  
  print STDERR "Alignment & phylogeny completed for: '$cluster'\n" unless $verbose;
  $pm->finish;	#end fork
}
$pm->wait_all_children;	


### subroutines
sub purge_names{
# purging fasta names of RAxML-unfriendly characters #
  my ($cluster, $curdir) = @_;
  
  my $cmd = "perl -pi -e \"s/[\\t :,)(\\]\\['#@]/_/g;s|/|_|g\"  $curdir/$prefix\_AA/$cluster";
  `$cmd`;
  
  $cmd = "perl -pi -e \"s/[\\t :,)(\\]\\['#@]/_/g;s|/|_|g\"  $curdir/$prefix\_nuc/$cluster";	
  `$cmd`;  
}

sub edit_stops{
  # changing any internal stops in AA sequences '*' to 'X'
  my ($cluster, $curdir) = @_;
  
  my $cmd = "perl -pi -e \"s/\\*/X/g\"  $curdir/$prefix\_AA/$cluster";
  `$cmd`;  
}



sub rm_PEGs{
  my ($cluster, $curdir, $prefix, $ML_dir, $rn_dir) = @_;
  
  my $tree_file = "$ML_dir/$cluster\_ML_lad.nwk";
  my $treeio = Bio::TreeIO -> new(-file => $tree_file,
								-format => 'newick');
  for my $tree ($treeio->next_tree){
    for my $node ($tree->get_nodes){
      next unless $node->is_Leaf;
      my $nodeID = $node->id;
      $nodeID =~ s/\.peg.+//g;
      $node->id($nodeID);
			}
    
    my $outfile = "$rn_dir/$cluster\_ML_lad_rn.nwk";
    my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
    $out->write_tree($tree);
    last;
  }
}

sub phy2raxml{
  my ($cluster, $curdir, $prefix, $pal2nal_dir, $ML_dir, $threads, $raxml_prog) = @_;	
	
  # convert to phylip #
  my $cmd = "alignIO.py $pal2nal_dir/$cluster $ML_dir/$cluster.phy";
  `$cmd`;
  
  # changing directories #
  chdir $ML_dir or die $!;
  
  # calling raxml #
  print STDERR "Starting RAxML inference on: '$cluster'\n";
  $cmd = "$raxml_prog -f a -x 0319 -p 0911 -# 100 -m GTRGAMMA -s $cluster.phy -n $cluster\_ML -T $threads";
  `$cmd`;
  
  if ($? == -1) {
    print "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  
  # ladderizing tree #
  rename("RAxML_bipartitions.$cluster\_ML", "$cluster\_ML.nwk") or die $!;
  $cmd = "ladderize.r -t $cluster\_ML.nwk";
  `$cmd`;
  
  if ($? == -1) {
    print "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  
  # moving back to current directory #
  chdir $curdir or die $!;
}

sub call_pal2nal{
  # calling mafft #
  my ($cluster, $curdir, $prefix, $align_dir, $pal2nal_dir, $rename_b, $pal2nal) = @_;
  die " ERROR: can't find nucleotide sequence!\n"
    unless -e  "$curdir/$prefix\_nuc/$cluster";
  `perl -pi -e 's/ /_/g' $align_dir/$cluster`;
  `perl -pi -e 's/ /_/g' $curdir/$prefix\_nuc/$cluster`;
  
  my $cmd = "$pal2nal $align_dir/$cluster $curdir/$prefix\_nuc/$cluster -output fasta";
  $cmd .= " | replaceOrgWithAbbrev.py" unless $rename_b;
  $cmd .= " > $pal2nal_dir/$cluster";
  `$cmd`;
  
  if ($? == -1) {
    print "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
}

sub call_mafft{
  # calling mafft #
  my ($cluster, $mafft_prog, $curdir, $prefix, $align_dir, $threads) = @_;
  
  # calling mafft #
  my $cmd = "$mafft_prog --quiet --thread $threads $curdir/$prefix\_AA/$cluster > $align_dir/$cluster";
  #print Dumper $cmd; exit;
  `$cmd`;
  
  if ($? == -1) {
    print "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    printf "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
}

sub get_seq_cnt{
  # count number of sequences in cluster #
  my ($cluster, $curdir) = @_; 
  
  open IN, "$curdir/$prefix\_AA/$cluster" or die $!;
  my $seq_cnt = 0;
  while(<IN>){
    $seq_cnt++ if /^>/;
  }
  close IN;
  
  return $seq_cnt;
}

sub get_fasta_files{
  my ($curdir, $prefix) = @_;
  
  opendir(DIR, "$curdir/$prefix\_AA");
  my @files = readdir(DIR);
  closedir DIR;
  
  @files = grep(! /^\./, @files);
  
  #print Dumper @files; exit;
  return \@files;
}

sub get_fastas{
  # loading list of clusters and making fasta files 
  my ($clusters_in, $curdir, $prefix) = @_;
  print STDERR "Getting AA fasta files from ITEP\n" unless $verbose;
  my $cmd_AA = "cat $clusters_in | db_getClusterGeneInformation.py | getClusterFastas.py $curdir/$prefix\_AA";
  `$cmd_AA`;
  print STDERR "Getting nuc fasta files from ITEP\n" unless $verbose;
  my $cmd_nuc = "cat $clusters_in | db_getClusterGeneInformation.py | getClusterFastas.py -n $curdir/$prefix\_nuc";
  `$cmd_nuc`;
  
  # checking for written fasta files
  opendir IN, "$curdir/$prefix\_AA" or die $!;
  die "ERROR: no AA files written! Is the right ITEP db sourced?\n" unless grep(!/^\./, readdir IN);
  closedir IN;
  
  opendir IN, "$curdir/$prefix\_nuc" or die $!;
  die "ERROR: no nuc files written! Is the right ITEP db sourced?\n" unless grep(!/^\./, readdir IN);
  closedir IN;
}

sub make_dirs{
  # making directories for:
  ## AA alignment
  ## pal2nal
  ## phy/raxmlk (rename in process)
  my ($curdir, $prefix) = @_;
  
  my $align_dir = "$curdir/$prefix\_AA_aln/";
  rmtree($align_dir) if -d $align_dir;
  mkdir $align_dir or die $!;
  
  my $pal2nal_dir = "$curdir/$prefix\_AA_aln_pal2nal/";
  rmtree($pal2nal_dir) if -d $pal2nal_dir;
  mkdir $pal2nal_dir or die $!;
  
  my $ML_dir = "$curdir/$prefix\_AA_aln_pal2nal_ML/";
  rmtree($ML_dir) if -d $ML_dir;
  mkdir $ML_dir or die $!;
  
  my $rn_dir = "$curdir/$prefix\_AA_aln_pal2nal_ML_rn/";
  rmtree($rn_dir) if -d $rn_dir;
  mkdir $rn_dir or die $!;
  
  return ($align_dir, $pal2nal_dir, $ML_dir, $rn_dir);
}





__END__

=pod

=head1 NAME

ITEP_align-tree.pl -- alignment & phylogeny inference pipeline

=head1 SYNOPSIS

ITEP_align-tree.pl [options] -c cluster_list.txt

=head2 Flags

=over

=item -clusters  <char>

List of clusters (runID\tclusterID).

=item -prefix  <char>

Output directory prefix. ["clusters"]

=item -rename  <bool>

Call replaceOrgWithAbbrev.py? [TRUE]

=item -threads  <int>

Number of threads to run for mafft & raxml. [2]

=item -fork  <int>

Number of clusters to process in parallel. [1]

=item -raxml  <char>

Which raxml to use? [raxmlHPC-PTHREADS-SSE3]

=item -mafft  <char>

Which mafft to use? [mafft-linsi]

=item -pal2nal  <char>

Name of pal2nal executable. [pal2nal.pl]

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_align-tree.pl

=head1 DESCRIPTION

Simple pipeline for aligning (mafft) and
inferring the phylogenies (RAxML) of ITEP
gene clusters.

=head1 EXAMPLES

=head2 Basic usage:

ITEP_align-tree.pl -c cluster_list.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITOL_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

