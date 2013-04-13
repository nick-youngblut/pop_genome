#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $group_in, $tmp_dir);
GetOptions(
	   "group=s" => \$group_in,
	   "directory=s" => \$tmp_dir,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
#die " ERROR: provide a group file in Mothur format!\n" unless $group_in;

foreach my $infile (@ARGV){
	die " ERROR: $infile not found!\n" unless -e $infile;
	$infile = File::Spec->rel2abs($infile);
	}

### MAIN
my $outdir = make_outdir($tmp_dir);

my %summary;
foreach my $infile (@ARGV){
	# converting fasta to table #
	my $tbl_file = fasta2txt($infile, $outdir);
	
	# running SNAP #
	call_SNAP($tbl_file);
	
	# parsing output #
	parse_SNAP_summary($tbl_file, $infile, \%summary);
	}

write_summary_table(\%summary);



### Subroutines
sub write_summary_table{
# writing summary table #
	my ($summary_r) = @_;
	
	foreach my $file (keys %$summary_r){
		print join("\t", $file, 
			$summary_r->{$file}{"dn"},
			$summary_r->{$file}{"ds"},
			$summary_r->{$file}{"dn/ds"}), "\n";
		}
	}

sub parse_SNAP_summary{
# parsing the output from SNAP #
	my ($tbl_file, $infile, $summary_r) = @_;
	
	open IN, "$tbl_file.summary" or die $!;
	while(<IN>){
		if (/Averages of all pairwise comparisons/){
			chomp;
			s/,//g;
			my @line = split / +/;	#7=ds, 10=dn, 13=dn/ds
			$summary_r->{$infile}{"ds"} = $line[7];
			$summary_r->{$infile}{"dn"} = $line[10];
			$summary_r->{$infile}{"dn/ds"} = $line[13];
			last;
			}
		}
	close IN;
	
	unless(exists $summary_r->{$infile}){
		$summary_r->{$infile}{"ds"} = "NA";
		$summary_r->{$infile}{"dn"} = "NA";
		$summary_r->{$infile}{"dn/ds"} = "NA";
		}
	}

sub call_SNAP{
# calling SNAP #
	my ($tbl_file) = @_;
	
	my $cmd;
	if($verbose){ $cmd = "SNAP.pl $tbl_file"; }
	else{ $cmd = "SNAP.pl $tbl_file 2>/dev/null"; }

	print STDERR "$cmd\n" if $verbose;
	`$cmd`;	
	}

sub fasta2txt{
# converting fasta to txt file #
	my ($infile, $outdir) = @_;
	open IN, $infile or die $!;

	# loading fasta #
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		$_ =~ s/#.+//;
 		next if  $_ =~ /^\s*$/;	
 		
 		if($_ =~ />.+/){
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{	
 			s/\./-/g;
 			$fasta{$tmpkey} .= $_; 
 			}
		}
	close IN;

	# writing fasta as table #
	my @parts = File::Spec->splitpath($infile);
	my $outfile = "$outdir/$parts[2]";
	$outfile =~ s/\.[^\.]+$|$/.txt/;
	
	open OUT, ">$outfile" or die $!;
	foreach my $name (keys %fasta){
		print OUT join("\t", $name, $fasta{$name}), "\n";
		}
	close OUT;
	
	return $outfile;
	}

sub make_outdir{
# making outdir for all of the temp files #
	my $outdir = shift;
	
	$outdir = "SNAP_batch_tmp" unless $outdir;
	rmtree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	
	$outdir = File::Spec->rel2abs($outdir);
	chdir $outdir or die $!;
	
	return $outdir;
	}


__END__

=pod

=head1 NAME

SNAP_batch.pl -- run SNAP.pl on many alignments

=head1 SYNOPSIS

SNAP_batch.pl [options] alignment(s).fna > SNAP_summary.txt

=head2 options

=over

=item -directory

Temporary file directory. ["SNAP_batch_tmp"]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc SNAP_batch.pl

=head1 DESCRIPTION

Run SNAP.pl on many gene alignments and get a summary of the output.

=head1 EXAMPLES

=head2 Usage:

SNAP_batch.pl [options] alignment(s).fna > SNAP_summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

