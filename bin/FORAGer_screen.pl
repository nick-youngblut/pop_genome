#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;
use List::Util qw/sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @clust_dirs, @contig_dirs);
my $fork = 0;
my $len_cutoff = 1;
GetOptions(
	   "clusters=s{,}" => \@clust_dirs,
	   "contigs=s{,}" => \@contig_dirs,
	   "length=f" => \$len_cutoff, 		# cutoff length for a contig (%)
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

	# blasting, filtering, writing results #
	foreach my $clust_file (keys %$files_r){
		$pm->start and next;
	
		# blasting #
		## tblastn cluster vs contig ##
		my $tblastn_r = tblastn_cluster_contig($clust_file, $files_r->{$clust_file}, $clust_dir, $contig_dir);
	
		## blastp cluster vs cluster ##
		my $blastp_min_bit = blastp_cluster($clust_file, $clust_dir);
	
		# filtering #
		## get stdev of cluster peg lengths ##
		my $clust_stdev = get_len_stdev($clust_file, $clust_dir);
	
		## filtering by score and length ##
		filter_blast($tblastn_r, $blastp_min_bit, $clust_stdev);
	
		$pm->finish;
		}
	}
$pm->wait_all_children;


### Subroutines
sub filter_blast{
# 
	}

sub get_len_stdev{
	my ($clust_file, $clust_dir) = @_;
	
	# getting lengths #
	open IN, "$clust_dir/$clust_file" or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if /^\s*$/;	
 		if(/^>/){
 			$fasta{$_} = 0;
 			$tmpkey = $_;	# changing key
 			}
 		else{  $fasta{$tmpkey} += length $_; }
		}
	close IN;
		#print Dumper %fasta; exit;

	# getting stdev of lengths #
	my $N = scalar keys %fasta;
	if($N > 1){ return stdev([values %fasta]); }		# if >1 peg in cluster
	else{ return (values %fasta)[0]; }
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

sub blast_norm{
	my ($tblastn_r, $blastp_r) = @_;
	
	#foreach my $q (keys $tblastn_r){
	#	
	#	}

	}

sub blastp_cluster{
# blastp of cluster vs cluster #
	my ($clust_file, $clust_dir) = @_;
	
	my $cmd = "blastp -subject $clust_dir/$clust_file -query $clust_dir/$clust_file -soft_masking true -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	print STDERR $cmd, "\n";
	open PIPE, "$cmd |" or die $!;	
	
	my $min_bit;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		$line[5] =~ s/\s+//g;
		if($min_bit){ $min_bit = $line[5] if $min_bit > $line[5]; }
		else{ $min_bit = $line[5]; }
		}
	close PIPE;

		#print Dumper $min_bit; exit;
	return $min_bit;	#query=>subject=>blast_res
	}

sub tblastn_cluster_contig{
# tblastn of contig vs cluster; parsing results #
	my ($clust_file, $contig_file, $clust_dir, $contig_dir) = @_;
	
	my $cmd = "tblastn -subject $contig_dir/$contig_file -query $clust_dir/$clust_file -soft_masking true -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	print STDERR $cmd, "\n";
	open PIPE, "$cmd |" or die $!;
	
	my %blast_res;
	while(<PIPE>){
		chomp;
			#print $_, "\n";
		my @line = split /\t/;
		$blast_res{$line[0]}{$line[1]} = [@line[2..$#line]];
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
			print STDERR " WARNING: no contig file found for $q!\n" unless @hits;
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

FORAGer_screen.pl [flags] 

=head2 Required flags

=item -contig

>=1 diretories of contig files produced by FORAGer_assemble.pl

=item -cluster

>=1 diretories of cluster files produced by FORAGer.pl

=head2 Optional flags

=over

=item -length

Length range cutoff (+/- the gene cluster length stdev * '-length). [1] 

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_screen.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head3 Length cutoff.   
Range is defined as the standard deviation of genes the target
cluster +/-(stdev_gene_cluster * '-length'). [1]

=head1 EXAMPLES

=head2 Usage method 1

FORAGer_screen.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

FORAGer_screen.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

