#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Set::IntervalTree;
use threads;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @sam_in, $database_file, $runID, $index_in);
my $gene_extend = 100;
my $clusterID_col = 2;
my $num_threads = 1;
my $outdir_name = "Mapped2GeneCluster";
GetOptions(
		"index=s" => \$index_in, 		# an index file of sam_file => FIG#
		"column=i" => \$clusterID_col,
		"runID=s" => \$runID,
		"extend=i" => \$gene_extend,	# bp to extend beyond gene (5' & 3')
		"outdir=s" => \$outdir_name, 	# name of output directory
		"threads=i" => \$num_threads,	
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a run ID!\n" unless $runID;
die " ERROR: provide an index file!\n" unless $index_in;

### MAIN
# loading info from ITEP #
my $clusterID_r = load_cluster_ids($clusterID_col);
my $gene_start_stop_r = load_gene_info($clusterID_r, $runID);

# pulling out reads mapping to each gene region #
my $index_r = load_index($index_in);

## threading ##
my @running;
my @threads;
my (%reads_mapped, %mapped_summary);
foreach my $sam_file (keys %$index_r){
	
	# checking for presence of genes in fig #
	unless (exists $gene_start_stop_r->{$index_r->{$sam_file}}){
		print STDERR " WARNING: no genes for FIG->", $index_r->{$sam_file}, ", skipping\n";
		next;
		}
	
	# finding mapped reads #
	my ($itrees_r, $reads_r) = load_interval_tree($sam_file);
	reads_mapped_to_region($sam_file, $index_r->{$sam_file}, 
			$gene_start_stop_r, $gene_extend, $itrees_r, $reads_r,
			\%reads_mapped, \%mapped_summary);	
	}

$outdir_name = make_outdir($outdir_name);
write_reads_mapped(\%reads_mapped, $outdir_name);
write_summary_table(\%mapped_summary, $outdir_name);

#---------- Subroutines -----------#
sub write_summary_table{
# writing summary table on reads mapped to FIGs for each gene cluster #
	my ($mapped_summary_r, $outdir_name) = @_;
	
	open OUT, ">$outdir_name/mapped_summary.txt" or die $!;
	
	# header #
	print OUT join("\t", qw/Cluster FIG N_reads Gene_length Reads_by_length/), "\n";
	
	# body #
	foreach my $cluster (keys %$mapped_summary_r){
		foreach my $fig (keys %{$mapped_summary_r->{$cluster}}){
			print OUT join("\t", $cluster, $fig, 
				$mapped_summary_r->{$cluster}{$fig}{"count"},
				$mapped_summary_r->{$cluster}{$fig}{"length"},
				$mapped_summary_r->{$cluster}{$fig}{"count"} /
				$mapped_summary_r->{$cluster}{$fig}{"length"}), "\n"; 
			}
		}
	
	close OUT;
	
	print STDERR "...Summary file written: $outdir_name/mapped_summary.txt\n";
	}

sub write_reads_mapped{
# writing out all reads mapped #
	my ($reads_mapped_r, $outdir_name) = @_;
	
	# status #
	print STDERR "...writing out read files to $outdir_name\n";
	
	foreach my $cluster (keys %$reads_mapped_r){
		(my $outfile = $cluster) =~ s/^/clust/;
		open OUTP, ">$outdir_name/$outfile\_FR.fna" or die $!;
		open OUTA, ">$outdir_name/$outfile\_A.fna" or die $!;
		
		foreach my $read (keys %{$reads_mapped_r->{$cluster}}){
			# paired-end reads #
			if(exists $reads_mapped_r->{$cluster}{$read}{1} && 
				exists $reads_mapped_r->{$cluster}{$read}{2}){ 
				foreach my $pair (sort keys %{$reads_mapped_r->{$cluster}{$read}}){
					foreach my $mapID (keys %{$reads_mapped_r->{$cluster}{$read}{$pair}}){
						print OUTP join("\n", ">$read $pair:", $reads_mapped_r->{$cluster}{$read}{$pair}{$mapID}), "\n";
						last; 		# just 1st map ID
						}
					}
				}
			# all reads #
			foreach my $pair (sort keys %{$reads_mapped_r->{$cluster}{$read}}){
				foreach my $mapID (keys %{$reads_mapped_r->{$cluster}{$read}{$pair}}){
					print OUTA join("\n", ">$read $pair:", $reads_mapped_r->{$cluster}{$read}{$pair}{$mapID}), "\n";
					last;			# just 1st mapID
					}
				}
			}
		close OUTP;
		close OUTA;
		}
	}

sub make_outdir{
# making output directory #
	my ($outdir_name) = @_;
	
	$outdir_name = File::Spec->rel2abs($outdir_name);

	rmtree($outdir_name) if -d $outdir_name;
	mkdir $outdir_name or die $!;
	
	return $outdir_name;
	}

sub reads_mapped_to_region{
# parsing out reads that mapped to each gene region #
## $gene_start_stop_r = fig=>cluster=>contig=>start/stop=>value
## foreach gene (start-stop): find reads mapped to gene (from itree) 
	my ($sam_file, $sam_fig, $gene_start_stop_r, 
		$gene_extend, $itrees_r, $reads_r, 
		$reads_mapped_r, $mapped_summary_r) = @_;

	my %warnings;
	foreach my $fig (keys %$gene_start_stop_r){			# FIG
		next unless $fig == $sam_fig;					# just looking for reads mapped to the same genome
		foreach my $cluster (keys %{$gene_start_stop_r->{$fig}}){		# gene cluster
			foreach my $contig (keys %{$gene_start_stop_r->{$fig}{$cluster}}){
				# contig of gene in gene tree? #
				unless(exists $itrees_r->{$contig}){
					print STDERR "WARNING: no reads mapped to FIG:$fig -> CONTIG:$contig\n"
						unless exists $warnings{$contig};		
					$warnings{$contig} = 1;				# so the warning only needs to be given once
					next;
					}
				
				my $res = $itrees_r->{$contig}->fetch(
						$gene_start_stop_r->{$fig}{$cluster}{$contig}{"start"} - $gene_extend,
						$gene_start_stop_r->{$fig}{$cluster}{$contig}{"stop"} + $gene_extend
						);
				
				unless(@$res){
					print STDERR "WARNING: no reads mapped to FIG:$fig -> Contig:$contig -> cluster:$cluster\n";
					next;
					}
				
				# loading read names #
				foreach my $id (@$res){		# read IDs
					$reads_mapped_r->{$cluster}{$$id[0]}{$$id[1]}{$$id[2]} = $$id[3];		# cluster->read_ID->pair->mapID = read
					}
					
				# loading summary (summing number of reads per gene) #
				$mapped_summary_r->{$cluster}{$fig}{"count"} += scalar @$res;
				$mapped_summary_r->{$cluster}{$fig}{"length"} =
					($gene_start_stop_r->{$fig}{$cluster}{$contig}{"stop"} + $gene_extend) - 
					($gene_start_stop_r->{$fig}{$cluster}{$contig}{"start"} - $gene_extend)
					unless exists $mapped_summary_r->{$cluster}{$fig}{"length"};
				}
			}
		}
		#print Dumper %$reads_mapped_r; exit;
	}

sub load_interval_tree{
	my ($sam_file) = @_;
	
	# status #
	print STDERR "...loading $sam_file\n";
	
	# loading reads as hash #
	open IN, $sam_file or die $!;
	my %reads;			# contig ->  read_name -> map_ID -> category -> value
	while(<IN>){
		chomp;
		next if /^@/; 	# skipping header
		next if /^\s*$/;	# skipping blank lines
		
		# parsing #
		my @line = split /\t/;
		## filtering ##
		next if $line[3] eq "" || $line[3] == 0;		# if not mapped
		
		## paired end info ##
		if($line[6] eq "="){ 		# if pair mapping
			my $line2 = <IN>;		# loading pair
			my @line2 = split /\t/, $line2;
			
			## check bitwise flag; read back to original orientation ##
			check_bitwise(\@line);
			check_bitwise(\@line2);
			
			## loading into hash ##
				# start & stop by + strand
			$reads{$line[2]}{$line[0]}{1}{$.}{"start"} = $line[3];						# contig->seq->start
			$reads{$line[2]}{$line[0]}{1}{$.}{"stop"} = $line[3] + length $line[9];		
			$reads{$line[2]}{$line[0]}{1}{$.}{"nuc"} = $line[9];	
			
			$reads{$line2[2]}{$line2[0]}{2}{$.}{"start"} = $line2[3];						# contig->seq->start
			$reads{$line2[2]}{$line2[0]}{2}{$.}{"stop"} = $line2[3] + length $line2[9];		
			$reads{$line2[2]}{$line2[0]}{2}{$.}{"nuc"} = $line2[9];			
			}	
		else{ 
			## check bitwise flag; read back to original orientation ##
			check_bitwise(\@line);		

			## loading into hash ##
			$reads{$line[2]}{$line[0]}{1}{$.}{"start"} = $line[3];						# contig->seq->start
			$reads{$line[2]}{$line[0]}{1}{$.}{"stop"} = $line[3] + length $line[9];		
			$reads{$line[2]}{$line[0]}{1}{$.}{"nuc"} = $line[9];	
			}
		}
	close IN;
	
	# loading interval tree #
	my %itrees;
	foreach my $contig (keys %reads){
		$itrees{$contig} = Set::IntervalTree->new;		# interval tree for each contig (position by contig)
		
		# inserting into interval trees #
		foreach my $read (keys %{$reads{$contig}}){
			foreach my $pair (keys %{$reads{$contig}{$read}}){
				foreach my $map_id (keys %{$reads{$contig}{$read}{$pair}}){
					$itrees{$contig}->insert(
							[$read, $pair, $map_id, $reads{$contig}{$read}{$pair}{$map_id}{"nuc"}],
							$reads{$contig}{$read}{$pair}{$map_id}{"start"} - 1, 
							$reads{$contig}{$read}{$pair}{$map_id}{"stop"} + 1);
					}
				}
			}
		}
	return \%itrees, \%reads;
	}
	
sub check_bitwise{
# checking bitwise flag in sam file line #
## rev or rev-comp sequence if needed ##
## returns read back to original orientation (for assembly) ##
## input = array-ref ##
	my ($arr_r) = @_;
	
	if($$arr_r[1] & 10){	# reverse
		$$arr_r[9] = reverse($$arr_r[9]);
		}
	elsif($$arr_r[1] & 20){	# rev-comp
		$$arr_r[9] = revcomp($$arr_r[9]);
		}
	}

sub revcomp{
        # reverse complements DNA #
        my $seq = shift;
        $seq = reverse($seq);
                #$seq =~ tr/[a-z]/[A-Z]/;
        $seq =~ tr/ACGTNBVDHKMRYSWacgtnbvdhkmrysw\.-/TGCANVBHDMKYRSWtgcanvbhdmkyrsw\.-/;
        return $seq;
        }

sub load_index{
# loading index file (sam => FIG#) #
	my ($index_in) = @_;
	
	open IN, $index_in or die $!;
	my %index;
	while(<IN>){
		chomp;
		s/#.+//;
		next if /^\s*$/;
		
		my @line = split /\t/;
		die " ERROR: $line[0] does not exist!\n" unless -e $line[0];
		$index{$line[0]} = $line[1];		# sam => fig
		}
	close IN;

		#print Dumper %index; exit;
	return \%index;
	}

sub load_gene_info{
# getting geneIDs for a cluster from ITEP #
	my ($clusterID_r, $runID) = @_;
	
	# status #
	print STDERR "...loading gene info from ITEP\n";
	
	# query ITEP #
	my @q;
	foreach my $clusterID (keys %$clusterID_r){
		push(@q, join("\\t", $runID, $clusterID));
		}
	my $q = join("\\n", @q);
	
	my $qq = "printf \"$q\" | db_getClusterGeneInformation.py |";
	
	# parsing ITEP output #
	open PIPE, $qq or die $!;
	my %gene_start_stop;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		
		# parsing #
		my @line = split /\t/;

		## loading start-stop ##
		(my $fig = $line[0]) =~ s/fig\||\.peg.+//g;			# FIG number
		$line[4] =~ s/^$line[2]\.//;
		if($line[7] eq "+"){
			$gene_start_stop{$fig}{$line[13]}{$line[4]}{"start"} = $line[5];
			$gene_start_stop{$fig}{$line[13]}{$line[4]}{"stop"} = $line[6];	
			}
		elsif($line[7] eq "-"){ 	# if neg strand, flipping
			$gene_start_stop{$fig}{$line[13]}{$line[4]}{"start"} = $line[6];
			$gene_start_stop{$fig}{$line[13]}{$line[4]}{"stop"} = $line[5];
			}
		else{ die " ERROR: 'strand' must be '+' or '-'\n"; }
		}
	
	close PIPE;
	
	# sanity check #
	die " ERROR: no gene information found for gene clusters!\n" if
		scalar keys %gene_start_stop == 0;
	
		#print Dumper %gene_start_stop; exit;
	return \%gene_start_stop;		# fig=>cluster=>contig=>start/stop=>value
	}

sub load_cluster_ids{
# loading cluster_id file in 'ITEP' format #
	my ($clusterID_col) = @_;
	
	my %clusterID;
	while(<>){
		chomp;
		s/#.+//;
		next if /^\s$/;
		
		my @line = split /\t/;
		$clusterID{$line[$clusterID_col -1]} = 1;
		}

		die " ERROR: no clusters provided!\n" if scalar keys %clusterID == 0;
		
		print STDERR "Number of clusters provided: ", scalar keys %clusterID, "\n"
			unless $verbose;

		#print Dumper %clusterID; exit;
	return \%clusterID;
	}


__END__

=pod

=head1 NAME

FindGeneClustersInReads.pl -- pulling out reads mapping to a gene cluster

=head1 SYNOPSIS

FindGeneClustersInReads.pl [flags] < clusterID_file.txt

=head2 Required flags

=over

=item -index

Index file listing *sam files & FIG IDs.

2 column (*txt): FILE.sam FIG

=item -runID

ITEP Run ID for the clusters of interest.

=back

=head2 Optional flags

=over

=item -outdir

Output directory name (location of all mapped read files). [./Mapped2GeneCluster/]

=item -extend

Number of base pairs to extend around the gene of interest (5' & 3'). [100]

=item -column

Cluster ID column in clusterID_file (starts at 1). [2]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FindGeneClustersInReads.pl

=head1 DESCRIPTION

Pull out all reads mapping to each gene (& region around the gene)
in each gene cluster.

=head2 output files

File prefix = 'clust#'

'_FR.fna' = Just paired-end reads that both mapped to the gene region

'_A.fna' = All reads that mapped to the gene region

=head2 SAM files

Reads from a genome of interest (e.g. a draft genome) should be mapped
onto all other genomes (producing the SAM files required).

You do not need to provide a *sam file for each FIG in each gene cluster;
however, genes from any genomes lacking *sam files will be skipped (ie. 
no reads pulled out for those genes; less total reads for the local assembly)

Multiple mappings are allowed for the same file (e.g. bowtie2 with '-k' flag).

=head2 Requires

ITEP must be setup so that it can be queried. 

=head1 EXAMPLES

=head2 Mapping reads to genomes

$ bowtie2-build genome1.fna genome1

$ bowtie2 -x genome1 -S genome1.sam -1 reads_F.fq -2 reads_R.fq -k 10

=head2 Make an index file

column1 = *sam file location

column2 = FIG ID

=head2 Get some gene clusters of interest

$ db_getAllClusterRuns.py | grep "mazei_I" | 
db_getClustersWithAnnotation.py "methyl coenzyme M reductase" >
mcr_clust.txt

=head2 Finding reads that mapped to each gene in the gene clusters of interest

FindGeneClustersInReads.pl -in index.txt -runID all_I_2.0_c_0.4_m_maxbit < mcr_clust.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

