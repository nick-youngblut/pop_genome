#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $qdir, $sdir, $regex, $list_in);
GetOptions(
	   "query=s" => \$qdir,			# query directory
	   "subject=s" => \$sdir,		# subject directory
	   "regex=s" => \$regex, 		# striping query name
	   "list=s" => \$list_in, 		# list of queries
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a query (read files) directory!\n" unless $qdir;
die " ERROR: provide a subject (reference genomes) directory!\n" unless $sdir;
map{ die " ERROR: $_ not found!\n" unless -d $_ } ($qdir, $sdir);
die " ERROR: provide a list of basenames (organism) for the query reads or a regex to strip off
the read pair info (e.g. '/_.+//')\n" unless $regex || $list_in;

$regex = load_regex($regex) if $regex;
die " ERROR: $list_in not found!\n" if $list_in && ! -e $list_in;

### MAIN
# loading directory files #
my $qfiles_r = load_dir($qdir);
my $sfiles_r = load_dir($sdir);

# loading list #
#my $org_list_r;
#if($list_in){ $org_list_r = load_list($list_in); }
#else{ $org_list_r = regex_sdir($regex, $qfiles_r); }
load_list($list_in)

# loading query files #
#$qfiles_r = get_read_files($org_list_r, $qfiles_r);


### Subroutines
sub get_read_files{
	my ($org_list_r, $qfiles_r) = @_;
	
	my %read_files;
	foreach my $org (@$ord_list_r){
		my @hits = grep(/$org/, @$qfiles_r);
		my $nhits = scalar @hits;
		die " ERROR: $org has $nhits read files! (need 1 or two)\n"
			unless scalar $nhits ==1 || $nhits ==2;
		
		}
	}	

sub regex_sdir{
	my ($regex, $qfiles_r) = @_;
	
	my %u;
	foreach my $qfile (@$qfiles_r){
		if($$regex[1]){ $qfile =~ s/$$regex[0]/$$regex[1]/; }
		else{ $qfile =~ s/$$regex[0]//; }
		$u{$qfile} = 1;
		}
	
		#print Dumper $qfiles_r; exit;
	return [keys %u];
	}

sub load_regex{
	my $regex = shift;
	$regex =~ s/^\///;
	$regex =~ s/\/$//;
	my @parts = split /\//, $regex;
	$parts[0] =~ qr/$parts[0]/;
	
	return \@parts;
	}

sub load_dir{
	my ($dir) = @_;
	opendir IN, $dir or die $!;
	my @files = grep(/\.(fq|fastq|fasta|fna|fa)$/, readdir IN);
	close IN;
	
	return \@files;
	}

sub load_list_old{
	my ($list_in) = @_;
	open IN, $list_in or die $!;

	my @org_list;
	while(<IN>){
		chomp;
		s/#.+//g;
		next if /^\s*$/;
		
		push(@org_list, $_);
		}
	close IN;
	
	return \@org_list;
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

