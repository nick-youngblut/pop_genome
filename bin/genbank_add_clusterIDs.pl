#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $organ_in, @genbank_in, $runID, $debug);
GetOptions(
	   "organisms=s" => \$organ_in,
	   "genbank=s{,}" => \@genbank_in,
	   "runID=s" => \$runID,
	   "debug" => \$debug,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

# I/O error & defaults #
die " ERROR: provide an organisms file!\n" unless $organ_in;
die " ERROR: provide >= genbank files!\n" unless @genbank_in;
map{ die " ERROR: $_ not found!\n" unless -e $_ } @genbank_in;
die " ERROR: provide a run ID!\n" unless $runID;

### workflow ###
# get fig from "organisms" file 
# echo "fig|1236903.88888.peg.17" | db_getClustersContainingGenes.py

#------#
# MAIN #
## loading organisms file ##
my $organ_r = load_organisms($organ_in);

## editing each genbank ##
foreach my $gen_in (@genbank_in){
	### status ###
	print STDERR "...processing $gen_in\n";
	
	### I/O ###
	my $seqo = Bio::SeqIO->new(-file => $gen_in, -format => "genbank")->next_seq;
	
	### getting fig ID ###
	my $figID = get_figID($seqo, $organ_r, $gen_in);

	### pulling clusterIDs from ITEP ###
	my $peg_clust_r;
	if($debug){		 # for debuging hits to ITEP db
		$peg_clust_r = get_clusterIDs_debug($seqo, $runID, $figID);
		}
	else{
		my ($tmpfile, $peg_cnt) = write_all_pegs($seqo, $figID);
		$peg_clust_r = get_clusterIDs($tmpfile, $runID, $peg_cnt);
		}
	
	### adding clusterIDs to genbank ###
	add_clusterIDs($peg_clust_r, $seqo, $runID);
	write_genbank($gen_in, $seqo);
	}

#-------------#
# Subroutines #
sub write_genbank{
	my ($gen_in, $seqo) = @_;

	(my $outfile = $gen_in) =~ s/\.[^\.]+$|$/_cID.gbk/;
	my $outio = Bio::SeqIO->new(-format => "genbank", -file => ">$outfile");
	
	$outio->write_seq($seqo);
	
	print STDERR "...Modified genbank written: $outfile\n";
	}

sub add_clusterIDs{
# adding clusterIDs to genbank file #
	my ($peg_clust_r, $seqo, $runID) = @_;
	
	my $peg_cnt = 0;
	for my $feat (grep {$_->primary_tag eq "CDS"} $seqo->get_SeqFeatures){
		$peg_cnt++;
		print STDERR " WARNING: peg$peg_cnt does not have a cluster ID!\n" unless exists $peg_clust_r->{$peg_cnt};
		
		$feat->add_tag_value("color", $peg_clust_r->{$peg_cnt}{$runID});
		}
	}

sub get_clusterIDs{
# getting clusterIDs for each peg #
	my ($tmpfile, $runID, $peg_cnt) = @_;
	
	my $cmd = "cat $tmpfile | db_getClustersContainingGenes.py | ";
	open IN, $cmd or die $!;
	my %peg_clust;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		next unless $line[0] eq $runID;
		$line[2] =~ s/.+peg\.//;
		$peg_clust{$line[2]}{$line[0]} = $line[1];		# fig.peg, run, cluster	
		}
	close IN;	
	
	# sanity check #
	my $Nhits = scalar keys %peg_clust;
	print STDERR " WARNING: number of pegs: $peg_cnt,  number of hits: $Nhits\n Are CRISPR spacers/repeats in the genbank file?\n"
		unless $Nhits == $peg_cnt;
		
	return \%peg_clust;
	}

sub write_all_pegs{
# writing all figs|pegs to file #
	my ($seqo, $figID) = @_;
	
	my $outfile = "runID_fig_peg.txt";
	open OUT, ">$outfile" or die $!;
	
	my $peg_cnt = 0;
	for my $feat (grep {$_->primary_tag eq "CDS"} $seqo->get_SeqFeatures){
		$peg_cnt++;
		print OUT "fig|$figID.peg.$peg_cnt\n";
		}
	close OUT;
	
	return $outfile, $peg_cnt;
	}

sub get_clusterIDs_debug{
# adding clusterIDs to genbank as '/color=' flag #
	my ($seqo, $runID, $figID) = @_;
	
	# foreach CDS #
	my %peg_clust;
	my $peg_cnt = 0;
	for my $feat (grep {$_->primary_tag eq "CDS"} $seqo->get_SeqFeatures){
		$peg_cnt++;
		my $tmpfile = write_pegs($figID, $peg_cnt);

		my $cmd = "cat $tmpfile | db_getClustersContainingGenes.py | ";
		open IN, $cmd or die $!;

		while(<IN>){
			chomp;
			my @line = split /\t/;
			next unless $line[0] eq $runID;
		$line[2] =~ s/.+peg\.//;
			$peg_clust{$line[2]}{$line[0]} = $line[1];		# fig.peg, run, cluster
			}
		close IN;		
		
		unless(exists $peg_clust{$peg_cnt}{$runID} ){
			print STDERR " WARNING 'fig|$figID.peg.$peg_cnt' had no hit to database for $runID\n";
			print Dumper $peg_clust{$peg_cnt};
			}
		}	
	 
	return \%peg_clust;
	}

sub write_pegs{
# writing out a file with fig|peg #
	my ($figID, $peg_cnt) = @_;
	
	my $outfile = "runID_fig_peg.txt";
	open OUT, ">$outfile" or die $!;
	my $cmd = "fig|$figID.peg.$peg_cnt\n";
	print STDERR $cmd;
	print OUT $cmd;
	close OUT;
	
	return $outfile;
	}

sub get_figID{
# getting figID from organism name #
	my ($seqo, $organ_r, $gen_in) = @_;

	foreach my $feat (grep {$_->primary_tag eq "source"} $seqo->get_SeqFeatures){
		my @org = $feat->get_tag_values("organism");
		die " ERROR: no organism tag in $gen_in\n" unless @org;
		die " ERROR: $org[0] not found in organism file!\n" unless exists $organ_r->{$org[0]};
		return $organ_r->{$org[0]};
		}
	}

sub load_organisms{
# loading organisms file #
	my ($organ_in) = @_;
	open IN, $organ_in or die $!;

	my %organ;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t+/;
		$organ{$line[0]} = $line[2];
		}
	close IN;
	
	return \%organ;
	}


__END__

=pod

=head1 NAME

genbank_add_clusterIDs.pl -- adding clusterIDs to each CDS (uses ITEP)

=head1 SYNOPSIS

genbank_add_clusterIDs.pl [flags]

=head2 required flags

=over

=item -organisms

ITEP 'organisms' file. Found in ITEP home directory.

=item -genbank

>=1 merged genbank file.

=item -runID

Clustering run ID. Use "db_getAllClusterRuns.py" to choose.

=back

=head2 optional flags

=over

=item -debug

Query ITEP one peg at a time to find pegs that do not return a hit 
to the run ID specified.

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc genbank_add_clusterIDs.pl

=head1 DESCRIPTION

Add cluster IDs to each CDS in 1 or more genbank files. 

Peg IDs are determined by ordering in the genbanks! The genbank files 
must be the same used for the ITEP database construction, except if 
they are not merged. If not merged, use EMBOSS 'union' to merge.

ITEP querying is done with 'db_getClustersContainingGenes.py'.

Cluster IDs are added to each CDS entry with the tag: '/color='

Output file names modifed from input "_cID.gbk".

=head1 EXAMPLES

=head2 Usage:

genbank_add_clusterIDs.pl -org organisms -g file_merged.gbk -r mazei_isolates_I_2.0_c_0.4_m_maxbit

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

