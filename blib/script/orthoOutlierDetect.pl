#!/usr/bin/env perl

=pod

=head1 NAME

orthoOutlierDetect.pl -- find potential outliers in gene clusters, possibly due to gene truncations

=head1 SYNOPSIS

orthoOutlierDetect.pl [flags] *fasta > outliers.txt

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -ceiling  <float>

Outlier must have lower seqID than '-ceiling' [99.9]

=item -floor  <int>

Outlier half must have > '-floor' gaps [4]

=item -x  <float>

Multiplier of stdev for calling outlier. [1]

=item -p  <int>

Number of files to process in parallel. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc orthoOutlierDetect.pl

=head1 DESCRIPTION

Find divergent sequences in an alignment,
which may be caused by a truncation or other
issue with the gene or alignment. These divergent
sequences may be false orthologs in a gene cluster.

Divergent sequences are identified as: 

(mean sequenceID of a sequence versus each other in alignment) < 
(mean sequenceID of the whole alignment - stdev sequenceID of whole
alignment * stdev_multiplier)

Possible truncations are identified by first spliting the alignment 
in half. Outliers in both seqID & gaps (similar to seqID outlier detection)
for 1 half but not the other will be flagged as possible truncations that
produced the observed sequence divergence in the alignment.

=head3 output columns

=over

=item * file_name

=item * outlier_sequence_name

=item * alignment_sequenceID_mean

=item * alignment_sequenceID_stdev

=item * outlier_sequenceID_mean

=item * truncation?

=item * truncation_half

=back

=head1 EXAMPLES

=head2 Basic usage:

orthoOutlierDetect.pl gene_cluster.fasta > outliers.txt

=head2 Many fasta files

orthoOutlierDetect.pl *.fasta > outliers.txt


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;
use Statistics::Descriptive;
use Parallel::ForkManager;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b);
my $ceiling = 99.9;
my $floor = 4;
my $stdev_x = 1;
my $fork = 0;
GetOptions(
	"ceiling=f" => \$ceiling,
	"floor=i" => \$floor,
	"x=f" => \$stdev_x,
	"processors=i" => \$fork,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error & defaults ---#
my @files = @ARGV;
map{ die "ERROR: cannot find file!\n" unless -e $_ } @files;


#--- MAIN ---#
my $pm = new Parallel::ForkManager($fork);
foreach my $file (@files){
	# forking #
	my $pid = $pm->start and next;


	# status #
	print STDERR "Processing: $file\n";	
	
	my $in  = Bio::AlignIO->new(-file   => $file,
                             -format => 'fasta');
    while(my $alno = $in->next_aln() ){
		my %output; 
		
	    # seqID outliers #
	    my $pIDs_r = pairwise_percent_identity($alno);                 
	    my ($mean,$stdev) = get_stats($pIDs_r);
	    my $means_r = get_mean_seqID_indiv($pIDs_r);
	    my $outliers_r = flag_outliers(\%output, $file, $means_r, $mean, $stdev, $stdev_x, $ceiling);
	    
	    # checking for a trucation #
	    if(%$outliers_r){		
	    	# split alignment #
	    	my @aln_halves = split_aln($alno);
			
			my %half_stats;
			for my $i (0..$#aln_halves){
				# status #
				print STDERR " Processing: $file -> half ", $i+1, "\n";
				
				my $aln_h = $aln_halves[$i];
				
				# gaps #
				my $gap_counts_r = gap_count($aln_h);                 
				my ($gap_mean, $gap_stdev) = get_gap_stats($gap_counts_r);
				flag_outliers_gaps(\%half_stats, $i, $gap_counts_r, $gap_mean, $gap_stdev, $stdev_x, $floor);
				
				
				# seqID #
				my $pIDs_r = pairwise_percent_identity($aln_h);                 
	    		my ($mean,$stdev) = get_stats($pIDs_r);
	    		my $means_r = get_mean_seqID_indiv($pIDs_r);
	 			flag_outliers_ret(\%half_stats, $i, $means_r, $mean, $stdev, $stdev_x, $ceiling);
				}
			flag_truncation(\%output, $file, \%half_stats);
	    	}
	    
	    write_output(\%output);
	    }
	
	$pm->finish;
	}
 $pm->wait_all_children;
 

#--- Subroutines ---#
sub write_output{
	my $output_r = shift;
	
	foreach my $id (keys %$output_r){
		foreach my $index (keys %{$output_r->{$id}}){
			map{ $_ = 'NA' unless defined $_ } @{$output_r->{$id}{$index}}[0..6];
			next if ${$output_r->{$id}{$index}}[0] eq "NA";
			print join("\t", @{$output_r->{$id}{$index}}), "\n";
			}
		}
	}

sub flag_truncation{
	my ($output_r, $file, $half_stats_r) = @_;
	
	foreach my $id (keys %$half_stats_r){
		foreach my $index (keys %{$half_stats_r->{$id}}){
			next unless exists $output_r->{$id}{$index};			# not an outlier sequence
			
			my $h1_out = 0;
			if(exists $half_stats_r->{$id}{$index}{0}){
				$h1_out = 1 if exists $half_stats_r->{$id}{$index}{0}{'seqID'}
							&& exists $half_stats_r->{$id}{$index}{0}{'gaps'};
				}
			my $h2_out = 0;
			if(exists $half_stats_r->{$id}{$index}{1}){
				$h2_out = 1 if exists $half_stats_r->{$id}{$index}{1}{'seqID'}
							&& exists $half_stats_r->{$id}{$index}{1}{'gaps'};			
				}
			
			if($h1_out == 1 && $h2_out == 0){
				push @{$output_r->{$id}{$index}}, ("Possible_truncation", 'first_half');
				}
			elsif($h2_out == 1 && $h1_out == 0){
				push @{$output_r->{$id}{$index}}, ("Possible_truncation", 'second_half');
				}
			else{
				push @{$output_r->{$id}{$index}}, ('NA', 'NA');
				}
			}
		}
	}

sub flag_outliers_gaps{
# finding sequences that are relatively divergent in the dataset #
## $mean = total mean
## $means_r = mean per seq
	my ($half_stats_r, $i, $gap_counts_r, $mean, $stdev, $stdev_x, $floor) = @_;

	#my %outliers;
	foreach my $id (keys %$gap_counts_r){
		foreach my $index (keys %{$gap_counts_r->{$id}}){
			next unless $gap_counts_r->{$id}{$index} > $floor;	# not outlier unless > floor
			if( $gap_counts_r->{$id}{$index} > $mean + $stdev * $stdev_x ){
				#$outliers{$id}{$index} = $gap_counts_r->{$id}{$index};
				$half_stats_r->{$id}{$index}{$i}{"gaps"} = $gap_counts_r->{$id}{$index};
				}
			}
		}
	#return \%outliers;
	}

sub get_gap_stats{
# getting mean & stdev of gap values #
	my ($gap_counts_r) = @_;
	
	my $stat = Statistics::Descriptive::Full->new();
	foreach my $id (keys %$gap_counts_r){
		foreach my $index (keys %{$gap_counts_r->{$id}}){
			$stat->add_data( $gap_counts_r->{$id}{$index} );
			}
		}
	return $stat->mean(), $stat->standard_deviation();
	}

sub gap_count{
# counting number of gaps per seq #
	my ($aln_r) = @_;
	
	my %gaps;
	my @seqs = $aln_r->each_seq(); 
	
	for my $i (0..$#seqs){
		my $gap_cnt=0;
		my $seq = $seqs[$i]->seq;
		$gap_cnt++ while ($seq =~ /-/g);
			#print Dumper $gap_cnt; exit;
		$gaps{$seqs[$i]->id}{$i} = $gap_cnt;
		}
		
		#print Dumper %gaps; exit;
	return \%gaps;
	}

sub flag_outliers_ret{
# finding sequences that are relatively divergent in the dataset #
## $mean = total mean
## $means_r = mean per seq
	my ($half_stats_r, $i, $means_r, $mean, $stdev, $stdev_x, $ceiling) = @_;

	#my %outliers;
	foreach my $id (keys %$means_r){
		foreach my $index (keys %{$means_r->{$id}}){
			next unless $means_r->{$id}{$index} < $ceiling;		# not outlier unless < ceiling
		
			if($means_r->{$id}{$index} < $mean - $stdev * $stdev_x ){
				#$outliers{$id}{$index} = $means_r->{$id}{$index};
				$half_stats_r->{$id}{$index}{$i}{"seqID"} = $means_r->{$id}{$index};
				}
			}
		}
	#return \%outliers;
	}

sub split_aln{
# spliting alignment into 2 halves #
	my ($alno) = @_;
	
	# getting alignment length & nseq #
	my ($nseq, $aln_len) = @_;
	foreach my $seq ( $alno->each_seq() ){
		$nseq++;
		$aln_len = $seq->length();
		}
	my $half = int($aln_len / 2);			# rounding down for middle of alignment (if odd)
	
	# splitting into halves
	my $alno_h1 = $alno->slice(1,$half);
	my $alno_h2 = $alno->slice($half+1,$aln_len);
	
		#print Dumper $alno_h1; exit;
	return $alno_h1, $alno_h2;
	}

sub flag_outliers{
# finding sequences that are relatively divergent in the dataset #
## $mean = total mean
## $means_r = mean per seq
	my ($output_r, $file, $means_r, $mean, $stdev, $stdev_x, $ceiling) = @_;

	my %outliers;
	foreach my $id (keys %$means_r){
		foreach my $index (keys %{$means_r->{$id}}){
			next unless $means_r->{$id}{$index} < $ceiling;		# not outlier unless < ceiling
		
			if($means_r->{$id}{$index} < $mean - $stdev * $stdev_x ){
				$output_r->{$id}{$index} = [$file, $id, 
						sprintf("%0.3f", $mean), 
						sprintf("%0.3f", $stdev * $stdev_x), 
						sprintf("%0.3f", $means_r->{$id}{$index})
						];
			
				$outliers{$id}{$index} = $means_r->{$id}{$index};
				}
			}
		}

	return \%outliers;
	}

sub get_mean_seqID_indiv{
# getting mean seq ID to 
	my ($pIDs_r) = @_;
	
	my %pIDs_indiv;
	foreach (@$pIDs_r){
		$pIDs_indiv{$$_[0]}{$$_[3]} = Statistics::Descriptive::Full->new()
			unless exists $pIDs_indiv{$$_[0]};
		$pIDs_indiv{$$_[0]}{$$_[3]}->add_data($$_[2]);
		
		$pIDs_indiv{$$_[1]}{$$_[4]} = Statistics::Descriptive::Full->new()
			unless exists $pIDs_indiv{$$_[1]};
		$pIDs_indiv{$$_[1]}{$$_[4]}->add_data($$_[2]);
		}
		
	my %means;
	foreach my $id (keys %pIDs_indiv){
		foreach my $index (keys %{$pIDs_indiv{$id}}){
			$means{$id}{$index} = $pIDs_indiv{$id}{$index}->mean();
			}
		}
		
		#print Dumper %means; exit;
	return \%means;
	}

sub get_stats{
	my ($pIDs_r) = @_;
	
	my $stat = Statistics::Descriptive::Full->new();
	foreach (@$pIDs_r){
		$stat->add_data( $$_[2] );
		}
	return $stat->mean(), $stat->standard_deviation();
	}

sub pairwise_percent_identity{
# getting array of pairwise identities # 
	my ($alno) = @_;
	
	#print Dumper $alno->length; exit;
	my @seqs = $alno->each_seq();
	
	my @pIDs;
	for my $i (0..$#seqs){
		for my $ii (0..$#seqs){
			next if $i >= $ii;
			my $alno_p = $alno->select($i+1, $ii+1);
			
			push @pIDs, [$seqs[$i]->id, $seqs[$ii]->id, $alno_p->percentage_identity, $i, $ii];
			#$pIDs{$seqs[$i]->id}{$seqs[$ii]->id}{"id"} = $alno_p->percentage_identity;
			#$pIDs{$seqs[$i]->id}{$seqs[$ii]->id}{"index"} = [$i, $ii];
			}
		}
	
		#print Dumper @pIDs; exit;
	return \@pIDs;
	}


