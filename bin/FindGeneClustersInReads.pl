#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Set::IntervalTree;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @sam_in, $database_file, $runID, $index_in);
my $gene_extend = 100;
my $clusterID_col = 2;
GetOptions(
		#"db=s" => \$database_file,		# database location
		#"sam=s{,}" => \@sam_in,			# all sam files for organisms of interest
		"index=s" => \$index_in, 		# an index file of sam_file => FIG#
		"column=i" => \$clusterID_col,
		"runID=s" => \$runID,
		"extend=i" => \$gene_extend,	# bp to extend beyond gene (5' & 3')
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a run ID!\n" unless $runID;
#die " ERROR: provide at least 1 sam file!\n" unless @sam_in;
#map{ die " ERROR $_ not found!\n" unless -e $_ } @sam_in;
die " ERROR: provide an index file!\n" unless $index_in;

### MAIN
### Connect to DB
#my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
#my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
#	or die " Can't connect to $database_file!";

# loading info from ITEP #
my $clusterID_r = load_cluster_ids($clusterID_col);
my $gene_start_stop_r = load_gene_info($clusterID_r, $runID);

	#foreach my $fig (keys %$gene_start_stop_r){				# FIG
	#	foreach my $cluster (keys %{$gene_start_stop_r->{$fig}}){	
	#		foreach my $contig (keys %{$gene_start_stop_r->{$fig}{$cluster}}){
	#			print $contig, "\n";
	#			}
	#		}
	#	}
	#exit;

# pulling out reads mapping to each gene region #
my $index_r = load_index($index_in);
foreach my $sam_file (keys %$index_r){
	# checking for presence of genes in fig #
	unless (exists $gene_start_stop_r->{$index_r->{$sam_file}}){
		print STDERR " WARNING: no genes for FIG->", $index_r->{$sam_file}, ", skipping\n";
		next;
		}
	
	# finding mapped reads #
	my ($itrees_r, $reads_r) = load_interval_tree($sam_file);
	reads_mapped_to_region($sam_file, $index_r->{$sam_file}, 
			$gene_start_stop_r, $gene_extend, $itrees_r, $reads_r);
			
	# 
	}


### Subroutines
sub reads_mapped_to_region{
# parsing out reads that mapped to each gene region #
## $gene_start_stop_r = fig=>cluster=>contig=>start/stop=>value
## foreach gene (start-stop): find reads mapped to gene (from itree) 
	my ($sam_file, $sam_fig, $gene_start_stop_r, 
		$gene_extend, $itrees_r, $reads_r) = @_;

	my %reads_mapped;		# reads mapped to a cluster region
	foreach my $fig (keys %$gene_start_stop_r){				# FIG
			#print Dumper "$fig -> $sam_fig";
		next unless $fig == $sam_fig;						# skipping if gene in other genome
		foreach my $cluster (keys %{$gene_start_stop_r->{$fig}}){		# gnee cluster
			foreach my $contig (keys %{$gene_start_stop_r->{$fig}{$cluster}}){
				# contig of gene in gene tree? #
				unless(exists $itrees_r->{$contig}){
					print STDERR "WARNING: no reads mapped to FIG:$fig -> CONTIG:$contig\n";
					next;
					}
				
				my $res = $itrees_r->{$contig}->fetch(
						$gene_start_stop_r->{$fig}{$cluster}{$contig}{"start"} - $gene_extend,
						$gene_start_stop_r->{$fig}{$cluster}{$contig}{"stop"} + $gene_extend
						);
				
				unless(@$res){
					print STDERR "\tWARNING: no reads mapped to  FIG:$fig -> Contig:$contig -> cluster:$cluster\n";
					next;
					}
				
				# loading read names #
				foreach my $id (@$res){
					#die " LOGIC ERROR: read does not exist: $id\n" unless
					#	exists $reads_r->{$id};
					#print Dumper $reads_r->{$contig}{$id};
					foreach my $read (keys %$id){
						$reads_mapped{$fig}{$cluster} = $res;	
						}
					}
				}
			}
		}
		exit;
		#print Dumper %reads_mapped; exit;
	return \%reads_mapped; 		# fig=>cluster=>[read_names]
	}

sub load_interval_tree{
	my ($sam_file) = @_;
	
	# status #
	print STDERR "...loading $sam_file\n";
	
	# loading reads as hash #
	open IN, $sam_file or die $!;
	my %reads;			# contig ->  seq -> ID# -> category -> value
	while(<IN>){
		chomp;
		next if /^@/; 	# skipping header
		next if /^\s*$/;	# skipping blank lines
		
		# parsing #
		my @line = split /\t/;
		## filtering ##
		next if $line[3] eq "" || $line[3] == 0;		# if not mapped
		
		## loading into hash ##
		$reads{$line[2]}{$.}{$line[0]}{"start"} = $line[3];						# contig->seq->start
		$reads{$line[2]}{$.}{$line[0]}{"stop"} = $line[3] + length $line[9];		
		$reads{$line[2]}{$.}{$line[0]}{"nuc"} = $line[9];
		
		### if paired read mapping ###		<-- skipping for now
		#if($line[8] =~ /\d+/){
		#	$reads{$line[2]}{$line[0]}{"pair_name"} = $line[9];
		#	}		
			
		last if $. > 1000000;
		}
	close IN;
	
	# loading interval tree #
	my %itrees;
	foreach my $contig (keys %reads){
		$itrees{$contig} = Set::IntervalTree->new;
		
		# inserting into interval trees #
		foreach my $read_id (keys %{$reads{$contig}}){
			foreach my $read (keys %{$reads{$contig}{$read_id}}){
				$itrees{$contig}->insert($read_id,
							$reads{$contig}{$read_id}{$read}{"start"} - 1, 
							$reads{$contig}{$read_id}{$read}{"stop"} + 1);
				}
			}
		}
	
	return \%itrees, \%reads;
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
		#print STDERR $q, "\n" unless $verbose;
	
	# parsing ITEP output #
	open PIPE, $qq or die $!;
	my %gene_start_stop;
	while(<PIPE>){
		chomp;
		next if /^\s*$/;
		
		# parsing #
		my @line = split /\t/;
			#print Dumper @line; exit;
			#print join("\t", $line[7], $line[10]), "\n"; 
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

		#$gene_info{$line[13]}{$line[0]}{$line[4]}{"nuc"} = $line[10];
		}
	
	close PIPE;
	
		#print Dumper %gene_start_stop; exit;
	return \%gene_start_stop;		# fig=>cluster=>contig=>start/stop=>value
	}

sub load_contigs{
# loading 'contig' table into memory #
	my ($dbh) = @_;
	
	my $q = "SELECT contig_mod, seq FROM contigs;";
	my $contig_seq_r = $dbh->selectall_arrayref($q); 
	
	my %contigs;
	foreach (@$contig_seq_r){
		$contigs{$$_[0]} = $$_[1];
		}
	
		#print Dumper %contigs; exit;
	return \%contigs;		# contig_mod => contig_seq
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
		
		print STDERR "Number of clusters provided: ", scalar keys %clusterID, "\n"
			unless $verbose;
			
		#print Dumper %clusterID; exit;
	return \%clusterID;
	}



__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

