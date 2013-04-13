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

my ($verbose, $names_in, $paup_block);
GetOptions(
	   "names=s" => \$names_in, 
	   "paup" => \$paup_block,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $names_r = load_names($names_in) if $names_in;
pres_abs_2_nexus($names_r);
write_paup_block() if $paup_block;

### Subroutines
sub load_names{
# laoding names of taxa to keep in table #
	my ($names_in) = @_;
	my %names;
	
	open IN, $names_in or die $!;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		$names{$line[0]}++; 
		}
	close IN;
	
	map{ print STDERR " WARNING: duplicate names provided!\n"
		if $names{$_} > 1 } keys %names;
		
	return \%names;
	}

sub pres_abs_2_nexus{
# converting pres-abs to nexus format readable by paup #
# removing 1st 3 lines #
# adding nexus header & footer

	my ($names_r) = @_;
	
	# loading file #
	my @pres_abs;
	my $ntaxa = 0;
	my $nchar = 0;
	while(<>){
		next unless $. > 3;		# skipping 1st 3 lines
		chomp;
		my @line = split /\t/;
		die " ERROR: matrix must be >= 2 columns\n" unless scalar @line >= 2;
		
		# checking taxon names #
		if($names_r && ! exists $names_r->{$line[0]}){
			print STDERR " $line[0] not in names file. Skipping.\n";
			next;
			}
			
		# loading #
		push @pres_abs, \@line;
		$ntaxa++;
		
		# nchar #
		die " ERROR: number of characters is not the same for $line[0]\n"
			if $nchar && $nchar != (scalar @line) -1;
		$nchar = (scalar @line) - 1;
		}

	# writing file #
	print "#NEXUS\n";
	print "begin data;\n";
	print "\t", join(" ", "dimensions", "ntax=$ntaxa", "nchar=$nchar"), ";\n";
	print "\tmatrix\n";
	
	foreach my $line (@pres_abs){
		print "\t\t", join("\t", $$line[0], join("", @$line[1..$#$line])), "\n";
		}
	print "\t\t;\n";
	print "\tend;";
	
	}
	
sub write_paup_block{
# writing paup block in script for running MP inference #
	
	print <<HERE;

begin paup;
	[outgroup #;  [designate the sequence numbers in your list of sequences]]
	set criterion = parsimony;
	set maxtrees=1 increase=no outroot=mono;
	hsearch;
	set maxtrees=100 increase=auto crit=parsimony;
	hsearch;
	savetrees root=yes brlens=yes file=PA_MP.tre replace=yes;
	quit;
end;	
HERE

	}


__END__

=pod

=head1 NAME

PresAbs2Nexus.pl -- convert ITEP P/A gene cluster table to nexus matrix

=head1 SYNOPSIS

PresAbs2Nexus.pl [options] < pres_abs.txt

=head2 options

=over

=item -name

*txt file with names to keep in table. 1st column used.

=item -paup

Include basic paup block for running an MP inference?
[FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc PresAbs2Nexus.pl

=head1 DESCRIPTION

Convert a gene cluster present absence table created with ITEP
to a nexus file with a character state matrix. 

=head1 EXAMPLES

=head2 Pipeing in P/A table from ITEP

db_getPresenceAbsenceTable.py -b | egrep "mazei_I_2.0_c_0.4_m_maxbit|runid" | transposeFile.py | PresAbs2Nexus.pl> mazei_PA.nex

=head2 Using certain taxa

=head3 getting names

echo mazei_I_2.0_c_0.4_m_maxbit | db_getOrganismsInClusterRun.py | perl -p -e 's/[ .]/_/g' > names.txt

=head3 Running PresAbs2Nexus.pl

db_getPresenceAbsenceTable.py -b | egrep "mazei_I_2.0_c_0.4_m_maxbit|runid" | transposeFile.py | PresAbs2Nexus.pl -n names.txt > mazei_PA.nex

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

