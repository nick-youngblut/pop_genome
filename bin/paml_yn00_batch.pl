#!/usr/bin/env perl

=pod

=head1 NAME

paml_yn00_batch.pl -- run SNAP.pl on many alignments and parse output

=head1 SYNOPSIS

paml_yn00_batch.pl [flags]

=head2 Required flags

=over

=item -directory  <char>

Directory containing fasta files (files extesion of "fasta", "fna", or "fa").

=item -group  <char>

Group file in Mothur format (2 column: 1st=taxon, 2nd=group).

=back 

=head2 Optional flags

=over

=item -prefix  <char>

Prefix of output files. ["yn00"]

=item -fork  <int>

Number of parallel yn00 calls. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc paml_yn00_batch.pl

=head1 DESCRIPTION

Run yn00 script from PAML package 
on a set of nucleotide alignments in 
fasta format.

A population (group) file must be provided 
for calculating mean omega (dN/dS) values
between populations. Format is 2 column: "TAXON_NAME\tPOPULATION".

Both the Nei & Gojobori 1986 method (used by SNAP)
and the Yang & Nielsen 2000 method
(incorporates base/codon frequency biases) are calculated.

=head2 Output files

=head3 PREFIX_results.txt

Mean pairwise omega values between 
sequences provided. 

=head3 PREFIX_summary.txt

Mean pairwise omega values between 
populations (excluding 'NA' values).

=head1 EXAMPLES

=head2 Basic usage:

paml_yn00_batch.pl -g groups.txt -d ./alignments/

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Parallel::ForkManager;
use IPC::Cmd qw/can_run/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $group_in, $fasta_dir);
my $prefix = "yn00";
my $fork = 1;
GetOptions(
	   "group=s" => \$group_in,
	   "directory=s" => \$fasta_dir,
	   "prefix=s" => \$prefix,			# output file prefix
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
can_run('yn00') or die 'yn00 is not installed!';
die "ERROR: provide a group file!\n" unless defined $group_in;
die "ERROR: provide a directory containing fasta files!\n" unless defined $fasta_dir;
die "ERROR: $fasta_dir not found!\n" unless -d $fasta_dir;
$fasta_dir = File::Spec->rel2abs($fasta_dir);


### MAIN
my $curdir = File::Spec->rel2abs(File::Spec->curdir());
my $files_r = load_file_list($fasta_dir, ["fasta", "fna", "fa"]);

my $group_r = load_group($group_in) if $group_in;

my $file_cnt = 0;
my $pm = Parallel::ForkManager->new($fork);
retreive_from_child($pm);
my %retrieved_responses;

# calling paml on each file #
foreach my $infile (@$files_r){
	# status #
	$file_cnt++;
	print STDERR "Number of files completed: $file_cnt\n" if $file_cnt % 100 == 0;
	my $pid = $pm->start and next;
		
	# tmp output dir #
	use File::Temp qw/tempdir/;
	my $tmpdir = File::Temp->newdir(); 		# temp directory
	my $outdir= $tmpdir->dirname;
	chdir $outdir or die $!;
	
	# symlink fasta #
	symlink("$fasta_dir/$infile", $infile) or die $!;
	
	# converting fasta to old-school phylip #
	my $fasta_r = read_fasta($infile);
	my ($outfile, $index_r) = fasta2phylip($fasta_r, $infile);
	
	# write out a yn00.ctl file #
	my ($ctlfile, $yn00_out) = write_ctl_file($outfile);
	
	# running yn00 (PAML) #
	`yn00 $ctlfile`;
	
	# parsing output #
	my $res_r = parse_yn00($yn00_out, $index_r, $group_r, $infile);
	
	#write_result_table($res_r, $prefix);
	
	# back to current directory
	chdir $curdir or die $!;
			
	$pm->finish(0, $res_r);
	}
$pm->wait_all_children;

# writing summary tables #
write_result_table(\%retrieved_responses, $prefix);
write_summary_table(\%retrieved_responses, $prefix);

### Subroutines
sub retreive_from_child{
	my $pm = shift;
	
	$pm -> run_on_finish (
		sub {
			#print Dumper @_; exit;
			my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
			
			#$ident = $$ unless defined $ident;
			
			# see what the child sent us, if anything
			if (defined($data_structure_reference)) {  # test rather than assume child sent anything
				my $reftype = ref($data_structure_reference);

				# we can also collect retrieved data structures for processing after all children have exited
				$retrieved_responses{$pid} = $data_structure_reference;
				} 
			else {
				print STDERR qq|ident "$ident" did not send anything.\n\n|;
				}
			}
		);
	}

sub write_summary_table{
# mean omega by group #
# summing by infile #
	my ($res_r, $prefix) = @_;
	
	open OUT, ">$prefix\_summary.txt" or die $!;
	print OUT join("\t", qw/method file group1 group2 mean_omega/), "\n";
	
	foreach my $pid (keys %$res_r){
		foreach my $method (keys %{$res_r->{$pid}}){
			foreach my $infile (keys %{$res_r->{$pid}{$method}}){
				my %omega_sum;
				my %omega_cnt;
				foreach my $row (@{$res_r->{$pid}{$method}{$infile}}){
					if(exists $omega_sum{$$row[0]}{$$row[1]}){
						if($$row[8] eq "NA"){
							$omega_sum{$$row[0]}{$$row[1]} += 0;
							$omega_cnt{$$row[0]}{$$row[1]} += 0;					
							}
						else{
							$omega_sum{$$row[0]}{$$row[1]} += $$row[8];
							$omega_cnt{$$row[0]}{$$row[1]}++;
							}
						}
					else{
						if($$row[8] eq "NA"){
							$omega_sum{$$row[1]}{$$row[0]} += 0;
							$omega_cnt{$$row[1]}{$$row[0]} += 0;
							}
						else{
							$omega_sum{$$row[1]}{$$row[0]} += $$row[8];
							$omega_cnt{$$row[1]}{$$row[0]}++;						
							}
						}
					}
				foreach my $group1 (keys %omega_sum){
					foreach my $group2 (keys %{$omega_sum{$group1}}){
						$omega_sum{$group1}{$group2} = 0 
							unless exists $omega_sum{$group1}{$group2};
						$omega_cnt{$group1}{$group2} = 0
							unless exists $omega_cnt{$group1}{$group2};
						my $mean;
						if($omega_cnt{$group1}{$group2} == 0){
							$mean = 'NA';
							}
						else{ $mean = $omega_sum{$group1}{$group2} / $omega_cnt{$group1}{$group2}; }
						print OUT join("\t", $method, $infile,
							$group1, $group2, $mean), "\n";
						}
					}
				}
			}
		}
	close OUT;	
	}

sub write_result_table{
# writing a table of all of the results #
	my ($res_r, $prefix) = @_;
	
	open OUT, ">$prefix\_results.txt" or die $!;

	my @header = qw/method file group1 group2
		seq1 seq2 S N t kappa omega dN+-SE dS+-SE/;
	
	print OUT join("\t", @header), "\n";
	foreach my $pid (keys %$res_r){
		foreach my $method (keys %{$res_r->{$pid}}){
			foreach my $infile (keys %{$res_r->{$pid}{$method}}){
				foreach my $row (@{$res_r->{$pid}{$method}{$infile}}){
					print OUT join("\t", $method, $infile, @$row), "\n";
					}
				}
			}
		}
	close OUT;
	}

sub parse_yn00{
### parsing yn00 result file ###
# parsing Nei & Gojobori; Yang & Nielsen #
	my ($outfile, $index_r, $group_r, $infile) = @_;
	
	open IN, "$outfile" or die $!;
	my %res;
	while(<IN>){
		if($_ =~ /^Nei & Gojobori 1986. dN\/dS \(dN, dS\)/){	#parsing Nei & Gojobori #
			my @dist;
			while(<IN>){
				last if /Yang & Nielsen \(2000\) method/;
				next if /^\s*$/;
				next if /^\S+ \S+/;
				push @dist, [split / {2,}| *\([^\)]+\) *|[\r\n]/];
				}
			# pairwise comparison table #
			my $dist_r = dist2tbl(\@dist);
			
			# loading hash #
			foreach my $row (@$dist_r){
				# getting group info #
				die "ERROR: cannot find group for '$_' in $infile\n"
					unless exists $group_r->{$index_r->{$$row[0]}};
				die "ERROR: cannot find group for '$_' in $infile\n"
					unless exists $group_r->{$index_r->{$$row[1]}};
				
				my $group1 = $group_r->{$index_r->{$$row[0]}};
				my $group2 = $group_r->{$index_r->{$$row[1]}};			
				
				# '-1.0000' values to 'NA' #
				$$row[2] = "NA" if $$row[2] == -1;
				
				# group1 group2 seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- #
				push @{$res{"NG"}{$infile}}, [$group1, $group2, @$row[0..1],
									("NA") x 4, $$row[2], ("NA") x 3 ];
				}
			}
		if(/^\(B\) Yang/){		# parsing Yang & Nielsen 2000
			my @dist;
			while(<IN>){ last if /^seq\. seq\./; }
			while(<IN>){
				chomp;
				last if /^\(C\) LWL85/;
				next if /^\s*$/;
				s/^ +//;
				my @line = split / +/;

				# getting group info #
				die "ERROR: cannot find group for '$_' in $infile\n"
					unless exists $group_r->{$index_r->{$line[0]}};
				die "ERROR: cannot find group for '$_' in $infile\n"
					unless exists $group_r->{$index_r->{$line[1]}};
				
				my $group1 = $group_r->{$index_r->{$line[0]}};
				my $group2 = $group_r->{$index_r->{$line[1]}};
				
				# values >= 99 to "NA" #
				$line[6] = "NA" if $line[6] >= 99;
				
				# group1 group2 seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- #
				push @{$res{"YN"}{$infile}}, [$group1, $group2, @line];
				}
			}
		}
	close IN;	
	
		#print Dumper %res; exit;
	return \%res;
	}

sub dist2tbl{
# converting lower triangle distance matrix to 2d array of pairwise comparisons #
	my ($dist_r) = @_;
	
	# getting max size #
	my $max = $#{$$dist_r[$#$dist_r]};	 
	
	my @dist;
	for my $i (0..$max){
		for my $ii (1..$max){
			#next if $i >= $ii; 				# lower triangle
			if(defined $$dist_r[$i][$ii]){
				push @dist, [$$dist_r[$i][0], $$dist_r[$ii][0], $$dist_r[$i][$ii]];
				}
			}
		}

		#print Dumper @dist; exit;
	return \@dist;
	}

sub write_ctl_file{
	my ($outfile) = @_;
	(my $yn00_out = $outfile) =~ s/\.[^.]+$|$/_res.txt/;
	(my $ctlfile = $outfile) =~ s/\.[^.]+$|$/.ctl/;

	open OUT, ">$ctlfile" or die $!;

	print OUT "     seqfile = $outfile      * sequence data file name
      outfile = $yn00_out     * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*       ndata = 1

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
";
	close OUT;
	
	return $ctlfile, $yn00_out;
	}

sub fasta2phylip{
# converting fasta 2 phylip & writing #
# making index of sequence order in output #
	my ($fasta_r, $infile) = @_;
	
	# number & length of sequences #
	my $nseq = scalar keys %$fasta_r;
	my $seq_len;
	foreach (keys %$fasta_r){
		if($seq_len){
			die "ERROR: sequences are not the same length in '$infile'\n"
				unless $seq_len = length($fasta_r->{$_});
			}
		else{ $seq_len = length($fasta_r->{$_}); }
		}
	
	# I/O #
	(my $outfile = $infile) =~ s/\.[^.]+$|$/_phy.txt/;
	open OUT, ">$outfile" or die $!;
	print OUT "   $nseq   $seq_len\n";

	# writing out sequences #
	my %index;
	my $cnt = 0;
	foreach my $seq (keys %$fasta_r){
		$cnt++;
		print OUT join("\n", "$cnt          ", $fasta_r->{$seq}), "\n";
		$index{$cnt} = $seq;
		}
	close OUT;
	
	return $outfile, \%index;
	}

sub read_fasta{
	my ($fasta_in) = @_;
	if(defined $fasta_in){
		die " ERROR: cannot find $fasta_in!" unless -e $fasta_in || -l $fasta_in;
		open IN, $fasta_in or die $!;
		}
	else{ *IN = *STDIN; }

	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			s/>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
	return \%fasta;
	}

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

sub load_file_list{
# loading list of files #
	my ($fasta_dir, $ext_r) = @_;
	
	opendir IN, $fasta_dir or die $!;
	my @files = readdir IN;
	closedir IN or die $!;
		#print Dumper @files; exit;
	my @fasta_files;
	foreach my $ext (@$ext_r){
		push @fasta_files, grep /\.$ext$/, @files;
		}
	
	# status #
	print STDERR "Number of fasta files found: ", scalar @fasta_files, "\n";
	
		#print Dumper @fasta_files; exit;
	return \@fasta_files;
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
		#print Dumper %group; exit;
	return \%group;
	}
