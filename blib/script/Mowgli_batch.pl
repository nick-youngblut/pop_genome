#!/usr/bin/env perl

=pod

=head1 NAME

Mowgli_batch.pl -- batch call and parse mowgli & its output, respectively

=head1 SYNOPSIS

Mowgli_batch.pl [flags] < input > output

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

Bootstrap cutoff for performing NNI by Mowgli

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc Mowgli_batch.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2

=head1 EXAMPLES

=head2 Basic usage:

Mowgli_batch.pl < input > output

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

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

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($species_tree, $gene_dir);
my ($verbose_b, $root_b, $DR_in);
my ($dup_cost, $trans_cost, $loss_cost) = (2,3,1);
my $boot_cutoff = 80;
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


#--- MAIN ---#
# getting list of gene trees #
my $gene_files_r = get_gene_tree_list($gene_dir);

# rerooting gene trees (if needed) #
## making temporary directory ##
my $tmpdir = File::Temp->newdir(); 		# temp directory
my $dirname = $tmpdir->dirname;

## rooting ##
my %file_index;
if(! $root_b){		# re-rooting
	foreach my $gene_file (@$gene_files_r){
		all_tree_rootings($gene_file, $dirname, \%file_index);
		}
	}
else{
	map{ $file_index{$_}{1} = "$gene_dir/$_"; } @$gene_files_r;
	}

# calling Mowgli #
my %res_all;
foreach my $gene_tree (keys %file_index){	
	my %res;
	foreach my $rooting (keys %{$file_index{$gene_tree}}){
		my $gene_tree_rooted = $file_index{$gene_tree}{$rooting};
		
		# loading gene tree values #
		$res{$gene_tree_rooted}{"species_tree"} = $species_tree;
		$res{$gene_tree_rooted}{"rooting"} = $rooting;
		$res{$gene_tree_rooted}{"dup_cost"} = $dup_cost;
		$res{$gene_tree_rooted}{"trans_cost"} = $trans_cost;
		$res{$gene_tree_rooted}{"loss_cost"} = $loss_cost;
		
		my $outdir = "$dirname/Mowgli_tmp_output_dir";
		call_Mowgli($species_tree, 
				$gene_tree_rooted, 
				$DR_in,
				$outdir,
				$dup_cost,
				$trans_cost,
				$loss_cost,
				$boot_cutoff, \%res);
				
		parse_costs($gene_tree_rooted, \%res, $outdir);
		parse_statistics($gene_tree_rooted, \%res, $outdir);
		}
	$res_all{$gene_tree} = \%res;
	}

# writing output #
write_table(\%res_all);


#--- Subroutines ---#
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
	print join("\t", "species_tree", "gene_tree", "gene_tree_rooted", @tree_cat, @DR_cat), "\n";

	# body #
	foreach my $gene_tree (keys %$res_all_r){
		foreach my $gene_tree_rooted (keys %{$res_all_r->{$gene_tree}}){
			# tree values #
			my @line = @{$res_all_r->{$gene_tree}{$gene_tree_rooted}}{@tree_cat};
			
			# donor-receiver values #
			foreach my $DR (keys %{$res_all_r->{$gene_tree}{$gene_tree_rooted}{"DR"}}){
				push @line, @{$res_all_r->{$gene_tree}{$gene_tree_rooted}{"DR"}{$DR}}{@DR_cat};
				}
			#print Dumper @line; exit;
			print join("\t", $gene_tree, $gene_tree_rooted, @line), "\n";
			}
		}
	#print Dumper %header; exit;
	}

sub parse_statistics{
### Description 
# parsing the statistics file from Mowgli
	my ($gene_tree_rooted, $res_r, $outdir) = @_;
	
	open IN, "gtlEnv_core_test/statistics.mpr" or die $!;
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
				$l[2] = 0 unless $l[2];
				$l[3] = 0 unless $l[3];
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
	-n 1 -T $boot_cutoff -f $DR_in -o $outdir 2>&1 ";
	$cmd =~ s/[\t\n]+/ /g;
	
		#print Dumper $cmd; exit;
	
	my $out = `$cmd`;
	die "Mowgli ERROR: $out!\n" if $out;
	}

sub all_tree_rootings{
#-- Description --#
# getting all possible tree rootings; just rooting by leaves #
	my ($gene_file, $tmpdir, $file_index_r) = @_;

	my $treeo = tree_io($gene_file);
	my @leaves = $treeo->get_leaf_nodes;

	# status #
	print STDERR "Re-rooting $gene_file...\n" unless $verbose_b;
	print STDERR "\tNumber of possible leaf rootings: ", 
		scalar @leaves, "\n" unless $verbose_b;
	
	my $root_cnt = 0;
	foreach my $leaf (@leaves){
		$root_cnt++;
		my $tree = tree_io($gene_file);
		my $reroot_leaf = $tree->find_node(-id => $leaf->id);
	
		$tree->reroot($reroot_leaf);	
		my $outfile = tree_write($tree, "$tmpdir/$gene_file", $root_cnt);
		$file_index_r->{$gene_file}{$root_cnt} = $outfile;
		
		# status #
		#print STDERR "\t\tCompleted rooting: $root_cnt\n" unless $verbose_b;
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

