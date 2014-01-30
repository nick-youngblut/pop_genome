#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $params, $extra_params);
my $assembler = "velvet";		# default assembler
my $in_dir;
GetOptions(
	   "directory=s" => \$in_dir,				# input directory
	   "assembler=s" => \$assembler, 	# assembler used
	   "parameters=s" => \$params, 		# assembly params
	   "extra=s" => \$extra_params, 	# params appended to default
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a directory!\n" unless $in_dir;
die " ERROR: $in_dir not found!\n" unless -d $in_dir;
$in_dir = File::Spec->rel2abs($in_dir);
chdir $in_dir or die $!;

$assembler = check_assembler_IO($assembler) unless $assembler eq "velvet";
$params = load_params($assembler) unless $params;
$params .= " $extra_params" if $extra_params;

### MAIN
my ($pair_files_r, $all_files_r) = find_read_files();
my %res;
if($assembler eq "idba_ud"){
	call_idba_ud($pair_files_r, $params, \%res);
	}
elsif($assembler eq "velvet"){
	call_velvet($all_files_r, $params, \%res);	
	}
else{ die " LOGIC ERROR: $!\n"; }

# writing out results #
write_file_list(\%res, $assembler);



### Subroutines
sub write_file_list{
# writing out files #
	my ($res_r, $assembler) = @_;
	
	foreach my $clust (keys %$res_r){
		foreach my $cat (keys %{$res_r->{$clust}}){
			print join("\t", $assembler, $clust, $cat, $res_r->{$clust}{$cat}), "\n";
			}
		}
	}

sub call_velvet{
	my ($all_files_r, $params, $res_r) = @_;
	
	foreach my $file (@$all_files_r){
		# sanity check #
		die $! unless -e $file;
		
		# editting params #
		$file =~ s/\.fna//;
		(my $tmp_params = $params) =~ s/%/$file/g;		# adding base name to params
		my @params = split /;/, $tmp_params;
		
		# calling velveth #
		my $cmd = "velveth $params[0]";
		print STDERR "$cmd\n" unless $verbose;
		`$cmd`;
		
		# calling velvetg #
		$cmd = "velvetg $params[1]";
		print STDERR "$cmd\n" unless $verbose;
		`$cmd`;
		
		# renaming contig file #
		eval{ rename("contigs.fa", "$file\_velvet-contig.fna") };
		
		# getting gene cluster #
		(my $clust = $file) =~ s/clust|_.+//g;
				
		# checking output #
		opendir IN, "." or die $!;
		my @contig = grep(/^$file\_velvet-contig.fna$/, readdir IN);
		close IN;		

		# loading results into hash #
		$res_r->{$clust}{"contigs"} = File::Spec->rel2abs($contig[0]) if $contig[0];
		}
	}

sub call_idba_ud{
	my ($pair_files_r, $params, $res_r) = @_;
	
	foreach my $file (@$pair_files_r){
		# sanity check #
		die $! unless -e $file;
		
		# editting params #
		$file =~ s/\.fna//;
		(my $tmp_params = $params) =~ s/%/$file/g;		# adding base name to params
		
		# calling idba_ud #
		my $cmd = "idba_ud $tmp_params";
		print STDERR $cmd, "\n" unless $verbose;
		`$cmd`;

		# getting gene cluster #
		(my $clust = $file) =~ s/clust|_.+//g;
		
		# checking output #
		opendir IN, "$file/" or die $!;
		my @files = grep(/^(contig.fa|scaffold.fa)$/, readdir IN);
		my @contig = grep(/contig/, @files); 
		my @scaf = grep(/scaffold/, @files);
		close IN;
		
		# loading results into hash #
		$res_r->{$clust}{"contigs"} = File::Spec->rel2abs($contig[0]) if $contig[0];
		$res_r->{$clust}{"scaffold"} = File::Spec->rel2abs($scaf[0]) if $scaf[0];
		}
	}

sub find_read_files{
# finding the read files produced by FORAGer.pl #
	opendir IN, "." or die $!;
	my @files = grep(/\.fna$/, readdir IN);
	close IN;
	
	my @pair_files = grep(/\_FR\.fna/, @files);
	my @all_files = grep(/_A\.fna/, @files);
		
		#print Dumper @all_files; exit;
	return \@pair_files, \@all_files;
	}

sub load_params{
# loading the default parameters for the assembler used #
	my ($assembler) = @_;
	my $params;
	if($assembler eq "idba_ud"){
		$params = "-r %.fna -o %";
		}
	elsif($assembler eq "velvet"){
		$params = ". 61 -fasta -short %.fna; . -cov_cutoff auto -exp_cov auto -ins_length 300 -ins_length_sd 100  -min_contig_lgth 100";
		}
	else{ die " LOGIG ERROR: $!\n"; }
	
	return $params;
	}

sub check_assembler_IO{
# checking the assembly #
	my ($assembler) = @_;
	if($assembler =~ /^i/i){ $assembler = "idba_ud"; }
	elsif($assembler =~ /^v/i){ $assembler = "velvet"; }
	else{ die " ERROR: assembler '$assembler' not recognized!\n"; }

	return $assembler;
	}

__END__

=pod

=head1 NAME

FORAGer_assemble.pl -- batch assembly of read files from FORAGer.pl

=head1 SYNOPSIS

FORAGer_assemble.pl [flags] > assembly_file_list.txt

=head2 Required flags

=over

=item -directory

Directory with the mapped read files from FORAGer.pl

=back

=head2 Optional flags

=over

=item -assembler

The assembler used on each file of reads (velvet | idba_ud). [velvet]

=item -parameters

The parameters provided to the assembler (see DESCRIPTION for defaults & details).

=item -extra

Extra parameters appended to the default parameters. 

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_assemble.pl

=head1 DESCRIPTION

Batch assemby of the read files produced by FORAGer.pl.
Each read file should contig all of the reads mapped to a gene cluster
(or just mapped paired-end reads).

The output is a *txt file containing the file locations for each assembly.
The columns are: 'assembler' 'gene_clusterID' 'contig|scaffold' 'file'

Paired-end reads used for idba_ud; all reads used for velvet.

=head2 Assembler parameters

=head3 Description

=over

=item '%'

 = The base file name (file name - extension of '.fna').
It will be replaced with the name of each read file.

=item ';'

 = Separates the parameters used for velveth & velvetg.

=back

=head3 Default idba_ud parameters

"-r %.fna -o %"

=head3 Default velvet parameters

". 61 -fasta -short %.fna; . -cov_cutoff auto -exp_cov auto -ins_length 300 -ins_length_sd 100  -min_contig_lgth 100"

=head2 WARNING

The script relies on standard file names "*_FR.fna" and "*_A.fna"
in the input file directory!

=head1 EXAMPLES

=head2 Basic usage (velvet assemblies)

FORAGer_assemble.pl > assembly_file_list.txt 

=head2 Using idba_ud

FORAGer_assemble.pl -a idba_ud > assembly_file_list.txt 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

