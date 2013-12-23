#!/usr/bin/env perl

=pod

=head1 NAME

bootstrap_byPop.pl -- get bootstrap stats by population (& all taxa)

=head1 SYNOPSIS

bootstrap_byPop.pl [flags] > bootstrap_stats.txt

=head2 Required flags

=over

=item -table  <char>

Table of tree files & matching gene cluster IDs (2-column; tab-delimited; 'file\tclusterID').

=item -population  <char>

Table of taxa in tree files & matching population (2-column; tab-delimited; 'taxon\tpopulation')

=back

=head2 Optional flags

=over

=item -format  <char>

Tree file format ('newick' or 'nexus'). ['newick']

=item -bootstrap  <bool>

Internal node IDs are bootstrap values? [TRUE]

=item -v  <bool>

Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc bootstrap_byPop.pl

=head1 DESCRIPTION

Get stats on how bootstrap values are distributed 
within and between populations. Stats will
also be provided for the entire tree ('total' population).

Not all taxa in the population file need to be found in each tree file.

=head1 EXAMPLES

=head2 Basic usage:

bootstrap_byPop.pl -p pop.txt -t tree-file-cluster_list.txt > boot-stats.txt

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
use Statistics::Descriptive;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $move_id_b, $pop_in, $file_in);
my $format = "newick";
GetOptions(
	"population=s" => \$pop_in,
	"table=s" => \$file_in,
	"format=s" => \$format,
	"bootstrap" => \$move_id_b,			# internal node id = bootstrap? [TRUE]
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error & defaults ---#
check_tree_format($format);
die "ERROR: provide a population file (-p)!\n" unless defined $pop_in;
die "ERROR: provide a file table (-t)!\n" unless defined $file_in;
die "ERROR: cannot find $file_in" unless -e $file_in;


#--- MAIN ---#
# loading file table #
my $files_r = load_file_table($file_in);

# loading population file #
my $pops_r = load_pop($pop_in) if $pop_in;		# taxon => population

# getting bootstrap stats for each population (& total) #
my $boot_r = get_bootstrap_stats($files_r, $pops_r, $format);

# writing out table #
write_stats_table($boot_r, $files_r);


#--- Subroutines ---#
sub write_stats_table{
	my ($boot_r, $files_r) = @_;
	
	my @cat = qw/min q1 mean median q3 max N/;
	
	print join("\t", qw/file cluster population/, @cat), "\n";
	foreach my $file (keys %$boot_r){
		foreach my $pop (keys %{$boot_r->{$file}}){
			# NA unless value for each stat category #
			map{ $boot_r->{$file}{$pop}{$_} = "NA" unless
					exists $boot_r->{$file}{$pop}{$_} } @cat;
					
			print join("\t", $file, $files_r->{$file}, $pop, @{$boot_r->{$file}{$pop}}{@cat}), "\n";
			}
		}
	}

sub get_bootstrap_stats{
# loading bootstraps for each tree and calculating stats #
	my ($files_r, $pops_r, $format) = @_;
	
	my %boot;
	foreach my $file (keys %$files_r){
		
		# status #
		print STDERR "...processing: $file\n"
			if $verbose_b;
				
		# loading tree #
		my $treeo = tree_io($file, $format);
		$treeo->move_id_to_bootstrap unless $move_id_b;
		
		# summing & counting bootstrap by population #
		## getting LCA of each pair of taxa ##
		my @pop_taxa = keys %$pops_r;
		my $boot_not_found = 0;
		my %stats;
		for my $i (0..$#pop_taxa){
			my @nodes1 = $treeo->find_node(-id => $pop_taxa[$i]);
			unless(defined $nodes1[0]){
				warn "WARNING: $pop_taxa[$i] not found in $file!\n"
					if $verbose_b;
				next;
				}
			
			for my $ii (0..$#pop_taxa){
				next if $i >= $ii;			# lower triange
	
				my @nodes2 = $treeo->find_node(-id => $pop_taxa[$ii]);
				unless(defined $nodes2[0]){
					warn "WARNING: $pop_taxa[$ii] not found in $file!\n"
						if $verbose_b;
					next;
					}

				my $lca = $treeo->get_lca(-nodes => [$nodes1[0], $nodes2[0]] );
				die "ERROR: no LCA for ", $nodes1[0]->id, " - ", $nodes2[0]->id, "\n"
					unless defined $lca;
				
				# checking for bootstrap; boostrap values of 0 will be missing, must be added back #
				unless($lca->bootstrap){
					$lca->bootstrap(0);
					$boot_not_found++;										
					}
				
				# adding bootstrap to stats object #
				my $pop_compare = join("__", 
							sort{$a cmp $b} ($pops_r->{$pop_taxa[$i]}, $pops_r->{$pop_taxa[$ii]} ) );

				## initializing if needed ##
				$stats{$pop_compare} = Statistics::Descriptive::Full->new()
					unless exists $stats{$pop_compare};
				$stats{'total'} = Statistics::Descriptive::Full->new()
					unless exists $stats{'total'};				

				$stats{$pop_compare}->add_data($lca->bootstrap);
				$stats{'total'}->add_data($lca->bootstrap);
				}
			}
			
		# warnings #
		warn "WARNING: $boot_not_found nodes in $file did not have bootstrap values (each given value of 0)!\n"
			if $boot_not_found > 0 && $verbose_b;
			
		# loading stats #
		foreach my $pop (keys %stats){
			my ($min, $q1, $mean, $median, $q3, $max, $N) = ('NA') x 7;
			$min = $stats{$pop}->min() if defined $stats{$pop}->min();
			$q1 = ($stats{$pop}->percentile(25))[0] if defined $stats{$pop}->percentile(25);
			$mean = $stats{$pop}->mean() if defined $stats{$pop}->mean();
			$median = $stats{$pop}->median() if defined $stats{$pop}->median();
			$q3 = ($stats{$pop}->percentile(75))[0] if defined $stats{$pop}->percentile(75);
			$max = $stats{$pop}->max() if defined $stats{$pop}->max();
			$N = $stats{$pop}->count() if defined $stats{$pop}->count();
			
			
			$boot{$file}{$pop} = {
				min => $min,
				q1 => $q1,
				mean => $mean,
				median => $median,
				q3 => $q3,
				max => $max,
				N => $N
				};
			}
		}
		
		#print Dumper %boot; exit;
	return \%boot;
	}

sub load_file_table{
# loading file table
# 	2 columns: file\tcluster
	my ($file_in) = @_;
	
	open IN, $file_in or die $!;
	my %files;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		die "ERROR: '$_' not in 2-column format!\n"
			unless scalar @l == 2;
		$files{$l[0]} = $l[1];
		}
	close IN;
	
		#print %files;
	return \%files;
	}

sub load_pop{
# loading population table #
	my ($pop_in) = @_;
	open IN, $pop_in or die $!;
	my %pops;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: population table is not formatted correctly\n"
			unless scalar @line == 2;
		$pops{$line[0]} = $line[1];
		}
	close IN;
		#print Dumper %pops; exit;
	return \%pops;
	}
	
sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}
	
sub tree_io{
	# loading tree object: just 1st tree #
	my ($tree_in, $format) = @_;
	my $input = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeio = $input->next_tree;	
	return $treeio;
	}




