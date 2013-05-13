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
use List::Util qw/max min/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @clust_dirs, @contig_dirs);
my $fork = 0;
my $len_cutoff = 1;
my $bit_cutoff = 0.4;
GetOptions(
	   "clusters=s{,}" => \@clust_dirs,
	   "contigs=s{,}" => \@contig_dirs,
	   "length=f" => \$len_cutoff, 			# cutoff length for a contig (%)
	   "bitscore=f" => \$bit_cutoff, 		# normalized bitscore cutoff
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide >=1 directory containing gene cluster fasta files (1 directory per query organism)!\n"
	unless @clust_dirs;
die " ERROR: provide >=1 directory containing FORAGer contig files (1 directory per query organism)!\n"
	unless @contig_dirs;
map{ die " ERROR: $_ not found!\n" unless -d $_; $_=File::Spec->rel2abs($_)} (@clust_dirs, @contig_dirs);


### MAIN
my $pm = new Parallel::ForkManager($fork);
for my $i (0..$#clust_dirs){
	# loading files #
	my $clust_dir = $clust_dirs[$i];
	my $contig_dir = $contig_dirs[$i];
	my $files_r = get_file_names($clust_dir, $contig_dir);
	die " ERROR: no files found!\n" unless %$files_r;

	# outdir #
	my $outdir = make_outdir($contig_dir, "_passed");

	# blasting, filtering, writing results #
	foreach my $clust_file (keys %$files_r){
		$pm->start and next;
		
		# blasting #
		## tblastn cluster vs contig ##
		my $tblastn_r = tblastn_cluster_contig($clust_file, $files_r->{$clust_file}, $clust_dir, $contig_dir);
	
		## blastp cluster vs cluster ##
		my $blastp_res_r = self_blast($clust_file, $clust_dir, "blastp");
		
		## blastn of contig vs contig ##
		#my $blastn_res_r = self_blast($files_r->{$clust_file}, $contig_dir, "blastn");
	
		## normalizing bit score: bit(contig vs gene) / max-bit(self-contig | self-gene) ##
		norm_blast($tblastn_r, $blastp_res_r);
	
		# loading fasta of contigs & clusters #
		my $contigs_r = load_fasta($files_r->{$clust_file}, $contig_dir);
		my $clusters_r = load_fasta($clust_file, $clust_dir);
	
		# filtering #
		my %summary;
		## filtering by number of tblastn hits (each gene must hit the contig) ##
		filter_by_Nhits(scalar keys %$clusters_r, $contigs_r, $tblastn_r, \%summary);
		
		## filtering by length cutoff ##	
		my $clust_range_r = get_clust_len_range($clusters_r);
		filter_by_length($clust_range_r, $len_cutoff, $contigs_r, $clusters_r, $tblastn_r, \%summary);	
	
		## filtering by score and length ##
		filter_by_bitscore($contigs_r, $tblastn_r, $bit_cutoff, \%summary);
	
		## writing out cluster fasta & PA table ##
		write_PA_table($clust_file, \%summary);
		write_passed_contig_fasta(\%summary, $contigs_r, $files_r->{$clust_file}, $outdir);
	
		$pm->finish;
		}
	}
$pm->wait_all_children;


### Subroutines
sub make_outdir{
# making output directory #
	my ($outdir_name, $append) = @_;
	
	$outdir_name .= "$append";
	
	$outdir_name = File::Spec->rel2abs($outdir_name);

	rmtree($outdir_name) if -d $outdir_name;
	mkdir $outdir_name or die $!;
	
	return $outdir_name;
	}

sub write_passed_contig_fasta{
# writing out a fasta of cluster & contig #
	my ($summary_r, $contigs_r, $contig_file, $outdir) = @_;

	(my $outfile = $contig_file) =~ s/\.[^\.]{1,6}$|$/_pass.fasta/; 
	open OUT, ">$outdir/$outfile" or die $!;

	foreach my $contig (keys %$summary_r){
		if($summary_r->{$contig}{"PA"}){		# passed; need to write out
			die " ERROR: $contig not found in contig file!\n" 
				unless exists $contigs_r->{$contig};
			print OUT join("\n", ">$contig", $contigs_r->{$contig}), "\n";
			}
		}
	
	close OUT;
	}
	
sub write_PA_table{
# writing out PA table to STDOUT #
	my ($clust_file, $summary_r) = @_;

	my @stats = qw/PA N_tblastn_hits_cutoff N_tblastn_hits length_cutoff hit_length_range norm_bit_score min_bit_score/;

	# writing body #
	foreach my $contig (keys %$summary_r){
		print join("\t", $clust_file, $contig);
		map{exists $summary_r->{$contig}{$_} ? print "\t$summary_r->{$contig}{$_}" : print "\tNA" } @stats;
		print "\n";
		}
	}

sub filter_by_bitscore{
# filtering by normalized bitscore #
	my ($contigs_r, $tblastn_r, $bit_cutoff, $summary_r) = @_;
	
	foreach my $subject (keys %$contigs_r){		# each contig
		if($bit_cutoff >= 0 && exists $tblastn_r->{$subject}){
			my @norm_bits;
			foreach my $query (keys %{$tblastn_r->{$subject}}){		# getting norm bit scores
				push(@norm_bits, ${$tblastn_r->{$subject}{$query}}[3]);
				}
			# pass/fail #
			my $x = min(@norm_bits);
			if($x >= $bit_cutoff){
				$summary_r->{$subject}{"norm_bit_score"} = "PASSED";
				}
			else{
				$summary_r->{$subject}{"PA"} = 0;
				$summary_r->{$subject}{"norm_bit_score"} = "FAILED";
				}
			$summary_r->{$subject}{"min_bit_score"} = sprintf("%.3f", $x);
			}
		else{
			$summary_r->{$subject}{"norm_bit_score"} = "NA";
			$summary_r->{$subject}{"min_bit_score"} = "NA";
			$summary_r->{$subject}{"PA"} = 0 if $bit_cutoff >= 0; 		# don't fail if not filtering by bit score
			}
		
		# PRESENT if passed all filters (no '0' for 'PA') #
		$summary_r->{$subject}{"PA"} = 1 unless exists $summary_r->{$subject}{"PA"};
		}
	
		#print Dumper %$summary_r; 
	}

sub filter_by_length{
# filtering the contigs by length relative to cluster #
	my ($clust_range_r, $len_cutoff, $contigs_r, $clusters_r, $tblastn_r, $summary_r) = @_;

	foreach my $contig (keys %$contigs_r){ 	
		if($len_cutoff && exists $tblastn_r->{$contig}){			# if tblastn hit(s)
			# checking all tblastn hit lengths #
			my $N_passed = 0;		# default = fail
			my @hit_lens;
			foreach my $query (keys %{$tblastn_r->{$contig}}){		# 1 hit per gene in cluster must meet length requirement
				my $hit_len = abs(${$tblastn_r->{$contig}{$query}}[1] - ${$tblastn_r->{$contig}{$query}}[0]);		# length in AA 
				$N_passed++ if $hit_len >= $$clust_range_r[0] - abs($$clust_range_r[0] - $$clust_range_r[0] * $len_cutoff)	# expanding negative range
								&& $hit_len <= $$clust_range_r[1] * $len_cutoff;			# expanding positive range
				push(@hit_lens, $hit_len);
				}
			## hit length range ##
			my $hit_range = join(":", sprintf("%.0f", min(@hit_lens)), 
									sprintf("%.0f", max(@hit_lens)));
			
			# pass/fail length cutoff #
			if($N_passed == scalar keys %$clusters_r){		# all hits must pass
				$summary_r->{$contig}{"length_cutoff"} = "PASSED";
				}
			else{
				$summary_r->{$contig}{"PA"} = 0;
				$summary_r->{$contig}{"length_cutoff"} = "FAILED";
				}
			$summary_r->{$contig}{"hit_length_range"} = $hit_range;
			}
		else{					# no tblastn hits at all
			$summary_r->{$contig}{"length_cutoff"} = "NA";
			$summary_r->{$contig}{"hit_length_range"} = "NA";			
			$summary_r->{$contig}{"PA"} = 0 if $len_cutoff; 		# no tblastn hits, so not PASSING
			}
		}
		#print Dumper $summary_r; exit;
	}

sub filter_by_Nhits{
# contig must be hit (tblastn) by all genes in cluster #
	my ($N_genes, $contigs_r, $tblastn_r, $summary_r) = @_;
	
	foreach my $contig (keys %$contigs_r){
		if(exists $tblastn_r->{$contig}){			# if tblastn hit(s)
			my $N_hits = scalar keys %{$tblastn_r->{$contig}};
			if($N_hits < $N_genes){
				$summary_r->{$contig}{"PA"} = 0;
				$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "FAILED";		
				}
			else{
				$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "PASSED";		
				}
			$summary_r->{$contig}{"N_tblastn_hits"} = $N_hits;		
			}
		else{
			$summary_r->{$contig}{"PA"} = 0;
			$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "FAILED";				
			$summary_r->{$contig}{"N_tblastn_hits"} = 0;		
			}
		}
	}

sub get_clust_len_range{
	my ($clusters_r) = @_;
	
	# getting lengths #
	my @lens;
	map{ push(@lens, length($_) * 3) } values %$clusters_r;

	# getting stdev of lengths #
	my $N = scalar @lens;
	if($N > 1){  return [min @lens,  max @lens]  }		# if >1 peg in cluster, max - min
	else{ return [$lens[0] * 0.75, $lens[0] * 1.25]; }	# just length of the gene +/- 0.25%
	}

sub stdev{
        my($data) = @_;
        
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub load_fasta{
	# version: 2.0
	# usage: load_fasta($fasta_file_name); returns a hash of sequences
	my ($fasta_in, $dir) = @_;
	open IN, "$dir/$fasta_in" or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper(%fasta); exit;
	return \%fasta;
	}

sub norm_blast{
# normalizing all tblastn bit scores: bit(contig vs gene) / max-bit(self-contig | self-gene) #
	my ($tblastn_r, $blastp_res_r, $blastn_res_r) = @_;
	
	foreach my $subject (keys %$tblastn_r){
			#die " ERROR: cannot find self blastn hit for $subject!\n"
			#	unless exists $blastn_res_r->{$subject};
		foreach my $query (keys %{$tblastn_r->{$subject}}){
			# sanity check #
			die " ERROR: cannot find self blastp hit for $query!\n"
				unless exists $blastp_res_r->{$query};

			#print join("\t", "tblastn: ${$tblastn_r->{$subject}{$query}}[3]", 
			#			"blastn: $blastn_res_r->{$subject}", 
			#			"blastp: $blastp_res_r->{$query}"), "\n";

			# normalizing #
			#${$tblastn_r->{$subject}{$query}}[3] = ${$tblastn_r->{$subject}{$query}}[3] / 
			#	max($blastn_res_r->{$subject}, $blastp_res_r->{$query});
			${$tblastn_r->{$subject}{$query}}[3] = ${$tblastn_r->{$subject}{$query}}[3] / 
				$blastp_res_r->{$query};
			}
		}
		#print Dumper %$tblastn_r; 
	}

sub self_blast{
# blasting against self; finding self hit w/ highest bit score #
	my ($file, $dir, $blast_type) = @_;
	
	my $cmd = "$blast_type -subject $dir/$file -query $dir/$file -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	$cmd .= " -soft_masking true" if $blast_type eq "blastp";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, "$cmd |" or die $!;	
	
	my %blast_clust;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		next unless $line[0] eq $line[1];
		$line[5] =~ s/\s+//g;
		if(exists $blast_clust{$line[0]}){ 		# getting max bit score if multiple self hits
			$blast_clust{$line[0]} = $line[5] if $line[5] > $blast_clust{$line[0]};
			}
		else{ 			# getting bit score if 1st self hit
			$blast_clust{$line[0]} = $line[5];
			}
		}
	close PIPE;

		#print Dumper $min_bit; exit;
	return \%blast_clust;	#query=>subject=>blast_res
	}

sub tblastn_cluster_contig{
# tblastn of contig vs cluster; parsing results #
	my ($clust_file, $contig_file, $clust_dir, $contig_dir) = @_;
	
	my $cmd = "tblastn -query $clust_dir/$clust_file -subject $contig_dir/$contig_file -soft_masking true -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, "$cmd |" or die $!;
	
	my %blast_res;
	while(<PIPE>){
		chomp;
			#print $_, "\n";
		my @line = split /\t/;
		if( exists $blast_res{$line[1]}{$line[0]} ){			
			$blast_res{$line[1]}{$line[0]} = [@line[2..$#line]] 
				if $line[5] > ${$blast_res{$line[1]}{$line[0]}}[3];
			}
		else{ $blast_res{$line[1]}{$line[0]} = [@line[2..$#line]]; }
		}
	close PIPE;
		#print Dumper %blast_res; exit;
	return \%blast_res;		# query=>subject=>blast_res
	}

sub get_file_names{
	my ($clust_dir, $contig_dir) = @_;

	# getting cluster files #
	opendir DIR, $clust_dir or die $!;
	my @clust_files = grep(/clust\d+.fasta/, readdir DIR);
	closedir DIR;
	
	# getting contig files #
	opendir DIR, $contig_dir or die $!;
	my @contig_files = grep(/clust\d+_A_velvet-contig.fasta|clust\d+_FR_idba_ud-contig.fasta/, readdir DIR);
	closedir DIR;
	
	#print Dumper @contig_files; exit;
	
	# making hash #
	my %files;
	foreach my $cluster (@clust_files){
		(my $q = $cluster) =~ s/\.[^.]+$//;
		my @hits = grep(/^$q\_/, @contig_files);

		if(@hits){	# if >= 1 hit
			die " ERROR: multiple contig files found for $q\n" if scalar @hits > 1;
			$files{$cluster} = $hits[0];
			}
		else{		# if contig file not found 
			print STDERR " WARNING: no contig file found for $q!\n" unless @hits || $verbose;
			$files{$cluster} = "NA";
			}
		}
		
		#print Dumper %files; exit;
	return \%files;		# cluster_file => contig_file
	}


__END__

=pod

=head1 NAME

FORAGer_screen.pl -- screening contigs to determine homology to the target gene cluster

=head1 SYNOPSIS

FORAGer_screen.pl [flags] > pres-abs_summary.txt

=head2 Required flags

=item -contig

>=1 diretories of contig files produced by FORAGer_assemble.pl

=item -cluster

>=1 diretories of cluster files produced by FORAGer.pl

=head2 Optional flags

=over

=item -bitscore

Normalized bit score cutoff (negative value to skip filtering). [0.4]

=item -length

Length range cutoff (+/- the gene cluster length stdev * '-length). 
'-length 0' skips filtering. [1]

=item -fork

Number of parallel cluster comparisons to perform. [1]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_screen.pl

=head1 DESCRIPTION

Determine which contigs produced by FORAGer_assemble.pl actually
should be included in the target gene cluster. Filtering is based on
sequence length & homology.

=head2 Normalized bit score cutoff

Normalized bit score = (tblastn_bit_score gene vs contig) / (blastp_bit_score gene vs gene).
The contig must hit all genes in cluster with a normalized bit score >= the cutoff.
The cutoff should probably be the same as used for the original gene clustering.

=head2 Length cutoff

The length range is defined as the min-max of gene lengths (bp) in the cluster.
'-length' is a scaling factor for how much to expand or contract this range.
By default, min-max length range values are used for the cutoff.

=head2 STDOUT



=head1 EXAMPLES

=head2 Usage:

FORAGer_screen.pl -contig neg_contigs/ -clust neg_clusters/ > pres-abs_summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

