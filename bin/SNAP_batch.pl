#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Forks::Super;
#use Parallel::ForkManager;
use IPC::Cmd qw/run can_run/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $group_in, $tmp_dir);
my $prefix = "SNAP_batch";
my $fork = 1;
my $SNAP_exec = 'SNAP.pl';
GetOptions(
	   "group=s" => \$group_in,
	   "directory=s" => \$tmp_dir,
	   "prefix=s" => \$prefix,			# output file prefix
	   "fork=i" => \$fork,
	   "SNAP_exec=s" => \$SNAP_exec,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
#die " ERROR: provide a group file in Mothur format!\n" unless $group_in;
can_run($SNAP_exec) or die "ERROR: cannot find $SNAP_exec in \$PATH\n";
foreach my $infile (@ARGV){
  die " ERROR: $infile not found!\n" unless -e $infile;
  $infile = File::Spec->rel2abs($infile);
}

### MAIN
my $group_r = load_group($group_in) if $group_in;
my $curdir = File::Spec->rel2abs(File::Spec->curdir());
#my $outdir = make_outdir($tmp_dir);



my %summary;
my %by_group;
foreach my $infile (@ARGV){
  my $job = fork{
    max_proc => $fork,
      share => [\%summary, \%by_group],
	sub => sub{
	  # tmp output dir #
	  use File::Temp qw/tempdir/;
	  my $tmpdir = File::Temp->newdir(); 		# temp directory
	  my $outdir= $tmpdir->dirname;
	  chdir $outdir or die $!;
	  
	  # converting fasta to table #
	  my $tbl_file = fasta2txt($infile, $outdir);
	  
	  # running SNAP #
	  call_SNAP($SNAP_exec, $tbl_file);
	  
	  # parsing output #
	  parse_SNAP_summary($tbl_file, $infile, \%summary, \%by_group, $group_r);
	  
	  # back to current directory
	  chdir $curdir or die $!;
	}
      }
}
waitall;

# writing summary tables #
write_by_group_table(\%by_group, $group_r, $prefix) if $group_r;
write_summary_table(\%summary, $prefix);

# clean up #
rm_fork_dirs();

### Subroutines
sub rm_fork_dirs{
# removing any fork directories still present #
  my $curdir = File::Spec->curdir();
  opendir IN, $curdir or die $!;
  my @files = readdir IN;
  closedir IN or die $!;
  
  foreach( grep(/^\.fhfork\d+/, @files) ){
    rmtree($_) or warn "Couldn't delete '$_'";
  }
  
  print STDERR "All Forks::Super tmp directories deleted\n";
}


sub write_by_group_table{
  # writing by-group summary of dn/ds #
  my ($by_group_r, $group_r, $prefix) = @_;
  
  # getting unique groups #
  my %ugroup;
  map{ $ugroup{$_} = 1 } values %$group_r;
  my @ugroup = keys %ugroup;
  
  # out FH #
  open OUT, ">$prefix\_by-group.txt" or die $!;
  
  # writing table #
  foreach my $file (keys %$by_group_r){
    for my $i (0..$#ugroup){
      for my $ii (0..$#ugroup){
	next if $i > $ii;		# lower triangle; keeping same comparisons
	
	# if no values in comparison that are not "NA" #
	if(! exists $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]} ||
	   $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"N"} =~ /NA/){
	  
	  print OUT join("\t", $file, $ugroup[$i], $ugroup[$ii], qw/NA NA NA NA/), "\n";
	  next;
	}
	
	# dn/ds #
	my $dn_ds = calc_dn_ds($by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"ds/dn"} /
			       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"N"});				
	
	# writing line (Ave values) #
	# header = file, group1, group2, ave_ds, ave_dn, ave_ds/dn, ave_dn/ds #
	print OUT join("\t", $file, $ugroup[$i], $ugroup[$ii],
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"ds"} /	
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"N"},
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"dn"} /
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"N"},
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"ds/dn"} /
		       $by_group_r->{$file}{$ugroup[$i]}{$ugroup[$ii]}{"N"},
		       $dn_ds), "\n";
      }
    }
  }
  
  close OUT;
}

sub write_summary_table{
  # writing summary table #
  my ($summary_r, $prefix) = @_;
  
  # out FH #
  open OUT, ">$prefix\_summary.txt" or die $!;
  
  # header: dn, ds, ds/dn, dn/ds	
  foreach my $file (keys %$summary_r){		
    print OUT join("\t", $file, 
		   $summary_r->{$file}{"ds"},
		   $summary_r->{$file}{"dn"},
		   $summary_r->{$file}{"ds/dn"},
		   $summary_r->{$file}{"dn/ds"}
		  ), "\n";
  }
  close OUT;
}

sub parse_SNAP_summary{
  # parsing the output from SNAP #
  my ($tbl_file, $infile, $summary_r, $by_group_r, $group_r) = @_;
  
  #print Dumper $group_r; exit;
  # getting group combinations #
  my %ugroup;
  map{ $ugroup{$_} = 1 } values %$group_r;
  my @ugroup = sort keys %ugroup;
  my @group_comb;
  for my $i (0..$#ugroup){
    for my $ii (0..$#ugroup){
      next if $ii < $i;
      push @group_comb, [$ugroup[$i], $ugroup[$ii]];
    }
  }
  
  open IN, "$tbl_file.summary" or die $!;
  while(<IN>){
    chomp;
    
    ## parsing by group (summing by group comparison) ##
    if (/^\d+ +\d+/ && $group_r){		# seq-seq comparison
      my @line = split / +/;
      next if $line[12] eq "NA";						# if no ds/dn
      
      die " ERROR: $line[2] not found in group file!\n"
	unless exists $group_r->{$line[2]};
      die " ERROR: $line[3] not found in group file!\n"
	unless exists $group_r->{$line[3]};
      
      ## loading by both grouping directions ##
      # ds #
      $by_group_r->{$infile}{$group_r->{$line[2]}} 			
	{$group_r->{$line[3]}}
	  {"ds"} += $line[10];
      
      $by_group_r->{$infile}{$group_r->{$line[3]}} 		
	{$group_r->{$line[2]}}
	  {"ds"} += $line[10];
      
      # dn #
      $by_group_r->{$infile}{$group_r->{$line[2]}} 			
	{$group_r->{$line[3]}}
	  {"dn"} += $line[11];
      $by_group_r->{$infile}{$group_r->{$line[3]}} 			
	{$group_r->{$line[2]}}
	  {"dn"} += $line[11];
      
      # ds/dn #
      $by_group_r->{$infile}{$group_r->{$line[2]}} 
	{$group_r->{$line[3]}}
	  {"ds/dn"} += $line[12];
      $by_group_r->{$infile}{$group_r->{$line[3]}} 
	{$group_r->{$line[2]}}
	  {"ds/dn"} += $line[12];
      
      # N ds/dn that are not "NA" #
      $by_group_r->{$infile}{$group_r->{$line[2]}} 
	{$group_r->{$line[3]}}
	  {"N"} += 1;
      $by_group_r->{$infile}{$group_r->{$line[3]}} 
	{$group_r->{$line[2]}}
	  {"N"} += 1;
    }
    
    ## mean values ##
    if (/Averages of all pairwise comparisons/){
      s/,//g;
      my @line = split / +/;	#7=ds, 10=dn, 13=dn/ds
      $summary_r->{$infile}{"ds"} = $line[7];
      $summary_r->{$infile}{"dn"} = $line[10];
      $summary_r->{$infile}{"ds/dn"} = $line[13];
      $summary_r->{$infile}{"dn/ds"} = calc_dn_ds($line[13]);
      last;
    }
  }
  close IN;
  
  # if no group values; NAs produced #
  foreach my $comb (@group_comb){
    unless(exists $by_group_r->{$infile}{$$comb[0]}{$$comb[1]}){
      $by_group_r->{$infile}{$$comb[0]}{$$comb[1]}{"ds"} = "NA";
      $by_group_r->{$infile}{$$comb[0]}{$$comb[1]}{"dn"} = "NA";
      $by_group_r->{$infile}{$$comb[0]}{$$comb[1]}{"ds/dn"} = "NA";
      $by_group_r->{$infile}{$$comb[0]}{$$comb[1]}{"N"} = "NA";									
    }
  }
  
  # if no mean values #
  unless(exists $summary_r->{$infile}){
    $summary_r->{$infile}{"ds"} = "NA";
    $summary_r->{$infile}{"dn"} = "NA";
    $summary_r->{$infile}{"ds/dn"} = "NA";
    $summary_r->{$infile}{"dn/ds"} = "NA";
  }
  
  #print Dumper %$by_group_r; exit; 
}

sub calc_dn_ds{
  # get dn/ds from ds/dn #
  my $ds_dn = shift;
  
  # dn_ds #
  my $dn_ds;
  if($ds_dn ne "NA" && $ds_dn > 0){
    $dn_ds = 1/$ds_dn;
  }
  else{ $dn_ds = "NA"; }	
  
  return $dn_ds;
}

sub call_SNAP{
# calling SNAP #
  my ($SNAP_exec, $tbl_file) = @_;
  
  my $cmd = "$SNAP_exec $tbl_file";
  #if($verbose){ $cmd = "$SNAP_exec $tbl_file"; }
  #else{ $cmd = "$SNAP_exec $tbl_file 2>/dev/null"; }
  
  print STDERR "$cmd\n" if $verbose;
  my( $success, $error_message, $full_buf,
      $stdout_buf, $stderr_buf ) =
        run( command => $cmd, verbose => 0 );
  
  unless ($success){
    die "ERROR: System call failed: '$cmd'\n" . join("\n", @$full_buf) . "\n";
  }

}

sub load_group{
# loading group file #
  my ($group_in) = @_;
  
  open IN, $group_in or die $!;
  my %group;
  while(<IN>){
    chomp;
    my @line = split /\t/;
    die " ERROR: group file must have >= 2 columns!\n"
      unless scalar @line >= 2;
    
    $group{$line[0]} = $line[1]; 
  }
  close IN;
  return \%group;
}

sub fasta2txt{
# converting fasta to txt file #
  my ($infile, $outdir) = @_;
  open IN, $infile or die $!;
  
  # loading fasta #
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    $_ =~ s/#.+//;
    next if  $_ =~ /^\s*$/;	
    
    if($_ =~ />.+/){
      $fasta{$_} = "";
      s/^>//;			# removing '>'
      $tmpkey = $_;	# changing key
    }
    else{	
      s/\./-/g;
      $fasta{$tmpkey} .= $_; 
    }
  }
  close IN;

  # writing fasta as table #
  my @parts = File::Spec->splitpath($infile);
  my $outfile = "$outdir/$parts[2]";
  $outfile =~ s/\.[^\.]+$|$/.txt/;
  
  open OUT, ">$outfile" or die $!;
  foreach my $name (keys %fasta){
    print OUT join("\t", $name, $fasta{$name}), "\n";
  }
  close OUT;
  
  return $outfile;
}

sub make_outdir{
  # making outdir for all of the temp files #
  my $outdir = shift;
  
  $outdir = "SNAP_batch_tmp" unless $outdir;
  rmtree($outdir) if -d $outdir;
  mkdir $outdir or die $!;
  
  $outdir = File::Spec->rel2abs($outdir);
  chdir $outdir or die $!;
  
  return $outdir;
}


__END__

=pod

=head1 NAME

SNAP_batch.pl -- run SNAP.pl on many alignments and parse output

=head1 SYNOPSIS

SNAP_batch.pl [options] alignment(s).fna

=head2 options

=over

=item -group  <char>

Group file in Mothur format (2 column: 1st=taxon, 2nd=group).
This is required to summarize dn/ds between groups.

=item -prefix  <char>

Prefix of output files. ["SNAP_batch"]

=item -directory  <char>

Temporary file directory. ["SNAP_batch_tmp"]

=item -fork  <int>

Number of parallel SNAP calls. [1]

=item -SNAP_exec  <char>

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc SNAP_batch.pl

=head1 DESCRIPTION

Run SNAP.pl on many gene alignments and get a summary of the output.
All values are averages of comparisons that did not produce 'NA' values.

One or more fasta files (nucleotide) can be provided as input.

If a group file is provided, the group file names (1st column) must
match the fasta names.

=head2 WARNINGS

Multi-copy genes should produce and output, but the usefulness of 
such output is questionable.

=head2 Output files

=head3 PREFIX_summary.txt

A summary for all comparisons in the alignment (excluding 'NA' values).
The columns are: file, ds, dn, ds/dn, dn/ds

=head3 PREFIX_by-group.txt

Only written if a group file is provided. The columns are:
file, group1, group2, ds, dn, ds/dn, dn/ds

=head1 EXAMPLES

=head2 Just averages for all taxa:

SNAP_batch.pl alignment(s).fna

=head2 Average values between each group

SNAP_batch.pl -g group.txt alignment(s).fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

