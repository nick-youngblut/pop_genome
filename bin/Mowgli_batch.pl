#!/usr/bin/env perl

=pod

=head1 NAME

Mowgli_batch.pl -- batch call and parse mowgli & its output, respectively

=head1 SYNOPSIS

Mowgli_batch.pl [flags] > output

=head2 Required flags

=over

=item -species  <char>

Species tree (newick; Mowgli format).

=item -gene  <char>

Directory of gene trees (all newick; Mowgli format).

=item -x  <char>

Donor-receiver file (Mowgli format). 

=back

=head2 Optional flags

=over

=item -root  <bool>

Mowgli runs on all possible rootings of each gene tree? [TRUE]

=item -D  <int>

Duplication cost. [2]

=item -T  <int>

Transfer cost. [3]

=item -L  <int>

Loss cost. [1]

=item -bootstrap  <int>

Bootstrap cutoff for performing NNI by Mowgli. [80]

=item -fork  <int>

Number of parallel Mowgli calls. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc Mowgli_batch.pl

=head1 DESCRIPTION

Calling Mowgli on multiple gene trees.
By default, Mowgli is called on
all possible leaf rootings of each gene
tree.

A donor-receiver file must be provided,
which specifies the clades that Mowgli
will sum the transfers to/from.

=head2 ERROR messages

"ERROR: child died with signal 6" caused by Mowgli error: "Exception: NNI() ----> invalid NNI".

=head1 EXAMPLES

=head2 Basic usage:

Mowgli_batch.pl -species species.nwk -gene gene_trees/ -x donor_receiver.txt > Mowgli_summary.txt

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
use Bio::TreeIO;
use File::Temp qw/ tempfile tempdir /;
use Forks::Super;
use File::Path;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($species_tree, $gene_dir);
my ($verbose_b, $root_b, $DR_in);
my ($dup_cost, $trans_cost, $loss_cost) = (2,3,1);
my $boot_cutoff = 80;
my $fork = 1;
GetOptions(
	   "species=s" => \$species_tree,
	   "gene_dir=s" => \$gene_dir,
	   "x=s" => \$DR_in,  
	   "root" => \$root_b,
	   "boot=i" => \$boot_cutoff,
	   "Duplication=i" => \$dup_cost,
	   "Transfer=i" => \$trans_cost,
	   "Loss=i" => \$loss_cost,
	   "verbose" => \$verbose_b,
	   "fork=i" => \$fork,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide a species tree!\n" unless defined $species_tree;
die "ERROR: provide a gene tree directory!\n" unless defined $gene_dir;
die "ERROR: provide a donor-receiver file!\n" unless defined $DR_in;
die "ERROR: cannot find $DR_in!\n" unless -e $DR_in;
die "ERROR: cannot find $species_tree!\n" unless -e $species_tree;
die "ERROR: cannot find $gene_dir!\n" unless -d $gene_dir;
$gene_dir = File::Spec->rel2abs($gene_dir);

# warnings #
print STDERR "WARNING: gene trees must be binary (use tree_multi2di.r)!\n";
print STDERR "\n";


#--- MAIN ---#
# getting list of gene trees #
my $gene_files_r = get_gene_tree_list($gene_dir);

# rerooting gene trees (if needed) #
## making directory for temporary files ##
my $dirname = "Mowgli_batch_tmp";
rmtree($dirname) or die $! if -d $dirname;
mkdir $dirname or die $!;
print STDERR "Writing temporary files to $dirname\n\n";

## rooting ##
my %file_index;
if(! $root_b){		# re-rooting
	foreach my $gene_file (@$gene_files_r){
		my $job = fork {
			share => [\%file_index],
			sub => sub{
				all_tree_rootings($gene_dir, $gene_file, $dirname, \%file_index);
				}
			};
		}
	waitall;
	}
else{
	map{ $file_index{$_}{1} = "$gene_dir/$_"; } @$gene_files_r;
	}

# calling Mowgli #
my %res_all;
foreach my $gene_tree (sort keys %file_index){	
  foreach my $rooting (sort{$a<=>$b} keys %{$file_index{$gene_tree}}){ 
   		my $gene_tree_rooted = $file_index{$gene_tree}{$rooting};
   		my $job = fork {
			max_proc => $fork,
  			share => [ \%res_all ],
    		sub => sub{
    			call_Mowgli_forked($dirname,$gene_tree,$gene_tree_rooted,$rooting);
    			}
    		};
    	}
	}
waitall;

# writing output #
write_table(\%res_all);


#--- Subroutines ---#
sub call_Mowgli_forked{
	my ($dirname, $gene_tree,$gene_tree_rooted,$rooting) = @_;
	
	# loading gene tree values #
	$res_all{$gene_tree}{$gene_tree_rooted}{"species_tree"} = $species_tree;
	$res_all{$gene_tree}{$gene_tree_rooted}{"rooting"} = $rooting;
	$res_all{$gene_tree}{$gene_tree_rooted}{"dup_cost"} = $dup_cost;
	$res_all{$gene_tree}{$gene_tree_rooted}{"trans_cost"} = $trans_cost;
	$res_all{$gene_tree}{$gene_tree_rooted}{"loss_cost"} = $loss_cost;

	my $outdir = "$dirname/Mowgli_tmp_output_dir";
		#my $outdir = "Mowgli_tmp_output_dir";		# debug

	call_Mowgli($species_tree, 
			$gene_tree_rooted, 
			$DR_in,
			$outdir,
			$dup_cost,
			$trans_cost,
			$loss_cost,
			$boot_cutoff, $res_all{$gene_tree});
	
	parse_costs($gene_tree_rooted, $res_all{$gene_tree}, $outdir);
	parse_statistics($gene_tree_rooted, $res_all{$gene_tree}, $outdir);
	}

sub write_table{
# writing output #
	my ($res_all_r) = @_;
	
	# header #
	my @tree_cat = qw/
		  species_tree
		  rooting
          dup_cost
		  trans_cost
		  loss_cost
          final_Cost
          number_of_Climbs
          number_of_bad_NNI
          number_of_tried_NNI/;
	my @DR_cat = qw/
		  Donor
		  Receiver
		  Tran
		  TranLoss/;        
	print join("\t", "gene_tree", "gene_tree_rooted", @tree_cat, @DR_cat), "\n";

	# body #
	foreach my $gene_tree (keys %$res_all_r){
		foreach my $gene_tree_rooted (keys %{$res_all_r->{$gene_tree}}){
			# tree values #
			my @line = @{$res_all_r->{$gene_tree}{$gene_tree_rooted}}{@tree_cat};
			map{ $_ = "NA" unless defined $_} @line;
			
			# donor-receiver values #
			foreach my $DR (keys %{$res_all_r->{$gene_tree}{$gene_tree_rooted}{"DR"}}){
				#push @line, @{$res_all_r->{$gene_tree}{$gene_tree_rooted}{"DR"}{$DR}}{@DR_cat};
				my @tmp = @line; 
				push @tmp, @{$res_all_r->{$gene_tree}{$gene_tree_rooted}{"DR"}{$DR}}{@DR_cat};
				map{ $_ = "NA" unless defined $_} @tmp[0..12];
				print join("\t", $gene_tree, $gene_tree_rooted, @tmp), "\n";
				}
			}
		}
	#print Dumper %header; exit;
	}

sub parse_statistics{
### Description 
# parsing the statistics file from Mowgli
	my ($gene_tree_rooted, $res_r, $outdir) = @_;
	
	open IN, "$outdir/statistics.mpr" or die $!;
	while(<IN>){
		chomp;
		if(/ \w.+:\d+$/){
			s/^ +//;
			my @l = split / +:/;
			$l[0] =~ s/ /_/g;
			$res_r->{$gene_tree_rooted}{$l[0]} = $l[1];
			}
		
		}
	close IN;
	
	#print Dumper %$res_r; exit;
	}

sub parse_costs{
### Description 
# parsing the cost file from Mowgli
# 	parsing number of trans & transLoss
	my ($gene_tree_rooted, $res_r, $outdir) = @_;
	
	open IN, "$outdir/costs.mpr" or die $!;
	
	while(<IN>){
		chomp;
		if(/^#Donor\tReceiver/){
			while(<IN>){
				chomp;
				last if /^\s*$/;
				my @l = split /\t+/;
				$l[0] = "NA" unless defined $l[0];
				$l[1] = "NA" unless defined $l[1];				
				$l[2] = 0 unless defined $l[2];
				$l[3] = 0 unless defined $l[3];
				
				# Donor-receiver name issues #
				$l[0] = "NA" if $l[1] eq "NA";			# both NA if 1 NA
				$l[1] = "NA" if $l[0] eq "NA";
				
				# loading hash #
				$res_r->{$gene_tree_rooted}{"DR"}{$.}{"Donor"} = $l[0];
				$res_r->{$gene_tree_rooted}{"DR"}{$.}{"Receiver"} = $l[1];
				$res_r->{$gene_tree_rooted}{"DR"}{$.}{"Tran"} = $l[2];
				$res_r->{$gene_tree_rooted}{"DR"}{$.}{"TranLoss"} = $l[3];
				}
			}
		}
	close IN; 
	
		#print Dumper %$res_r; exit;
	}

sub call_Mowgli{
### Description
# calling Mowgli #
	my ($species_tree, $gene_tree, $DR_in, $outdir,
		$dup_cost, $trans_cost, $loss_cost, $boot_cutoff, $res_r) = @_;
	
	print STDERR "Calling Mowgli on: '$gene_tree'\n" unless $verbose_b;
	
	my $cmd = "Mowgli_linux_i386 -s $species_tree -g $gene_tree 
	-d $dup_cost -t $trans_cost -l $loss_cost 
	-n 1 -T $boot_cutoff -f $DR_in -o $outdir";
	$cmd =~ s/[\t\n]+/ /g;
	
	system($cmd);

	if ($? == -1) {
    	die "ERROR: failed to execute: $!\n";
		}
	elsif ($? & 127) {
    	printf STDERR "ERROR: child died with signal %d, %s coredump\n",
        	($? & 127),  ($? & 128) ? 'with' : 'without';
        exit(1);
		}
	elsif( $? != 0 ) {
    	printf STDERR "Mowgli exited with value %d\n", $? >> 8;
		}
	
	}

sub all_tree_rootings{
### Description ###
# getting all possible tree rootings; just rooting by leaves #
	my ($gene_dir, $gene_file, $tmpdir, $file_index_r) = @_;

	# getting leaves #
	my $treeo = tree_io("$gene_dir/$gene_file");
	my @leaves = $treeo->get_leaf_nodes;

	# status #
	print STDERR "Using nw_reroot to re-root $gene_file...\n" unless $verbose_b;
	print STDERR "\tNumber of possible leaf rootings: ", 
		scalar @leaves, "\n" unless $verbose_b;
	
	# calling nw_reroot to reroot trees #
	my $root_cnt = 0;
	foreach my $leaf (@leaves){
		$root_cnt++;
		#my $tree = tree_io("$gene_dir/$gene_file");
		(my $outfile = $gene_file) =~ s/\.[^.]+$|$/_r$root_cnt.nwk/;
		my $cmd = join(" ", "nw_reroot $gene_dir/$gene_file",
						$leaf->id, ">$tmpdir/$outfile");				
		
		system($cmd);
		
		if ($? == -1) {
    		die "ERROR: failed to execute: $!\n";
			}
		elsif ($? & 127) {
    		printf "ERROR: child died with signal %d, %s coredump\n",
        		($? & 127),  ($? & 128) ? 'with' : 'without';
	        exit(1);
			}
		elsif($? != 0) {
    		printf STDERR "nw_reroot exited with value %d\n", $? >> 8
			}
			
		# loading file index #
		$file_index_r->{$gene_file}{$root_cnt} = "$tmpdir/$outfile";
		}
		
		#print Dumper $file_index_r; exit;
	}
	
sub tree_write{
### writting out a newick tree file ###
	my ($treeo, $outfile, $root_cnt) = @_;
	$outfile =~ s/\.[^.]+$|$/_r$root_cnt.nwk/;
	my $out = new Bio::TreeIO(-file => ">$outfile", -format => "newick");
	$out->write_tree($treeo);
	return $outfile;
	}

sub get_gene_tree_list{
#-- Description --#
# getting a list of gene trees specified by user #
	my ($gene_dir) = @_;
	
	opendir IN, $gene_dir or die $!;
	my @files = grep(!/^\./, readdir IN);
	closedir IN;
	
	# sanity check #
	die "ERROR: no gene trees found in '$gene_dir'!\n"
		unless scalar @files >0;
		
		#print Dumper @files;  exit;
	return \@files;
	}
	
sub tree_io{
	# loading tree object #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeo = $input->next_tree;	
	return $treeo;
	}

