#!/usr/bin/env perl

=pod

=head1 NAME

ITEP_organismsByGroup.pl -- parse and ITEP organisms file by groups in the ITEP groups file

=head1 SYNOPSIS

ITEP_organismsByGroup.pl [flags] 

=head2 Required flags

=over

=item -org

organisms file

=item -group

groups file

=back

=head2 Optional flags

=over

=item -dir

Output directory [organisms_byGroup].

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc ITEP_organismsByGroup.pl

=head1 DESCRIPTION

Parse and ITEP organisms file by groups in the ITEP groups file.
Can be useful for having ITEP-independent way to identifying taxa
in each designate groups.

=head1 EXAMPLES

=head2 Basic usage:

ITEP_organismsByGroup.pl -o organisms -g group

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
use File::Path qw/rmtree/;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $org_in, $group_in);
my $outdir = "organisms_byGroup";
GetOptions(
	"directory=s" => \$outdir,
	"organisms=s" => \$org_in,
	"groups=s" => \$group_in,
	"verbose" => \$verbose_b,
	"help|?" => \&pod2usage # Help
	);

#--- I/O error ---#
check_file($org_in, 'organisms');
check_file($group_in, 'group');

#--- MAIN ---#
# loading input #
my $org_r = load_org($org_in);
my $group_r = load_group($group_in);

# making output dir #
make_outdir($outdir);

# parsing organism file and writing to individual files #
parse_org($org_r, $group_r, $outdir);


#--- Subroutines ---#
sub parse_org{
	my ($org_r, $group_r, $outdir) = @_;
	
	foreach my $group (keys %$group_r){
		# making output file #
		open OUT, ">$outdir/$group" or die $!;
		
		foreach my $taxon (@{$group_r->{$group}}){
			if(exists $org_r->{$taxon}){
				print OUT join("\t", $taxon, $org_r->{$taxon}), "\n";
				}
			else{
				warn "WARNING: group->$group, taxon->$taxon not found in organism file!\n";
				next;
				}
			}
		close OUT;
		}
	}

sub make_outdir{
	my ($outdir) = @_;
	# removing existing directory #
	rmtree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	}


sub load_group{
	my ($group_in) = @_;
	
	open IN, $group_in or die $!;
	my %group;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /\t|;/;
		die "ERROR: line $. in group file is not formated correctly!\n"
			unless scalar @l >= 2;
		$group{$l[0]} = [@l[1..$#l]];
		}
	close IN;
		
	return \%group;
	}

sub load_org{
	my ($org_in) = @_;
	open IN, $org_in or die $!;
	my %org;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /\t/;
		die "ERROR: organisms file < 2 columns (tab-delimited)\n"
			unless scalar @l >= 2;
		
		$org{$l[0]} = $l[1];
		}
	close IN;
	return \%org;
	}

sub check_file{
	my ($infile, $name) = @_;
	die "ERROR: no $name file provided!\n"
		unless defined $infile;
	die "ERROR cannot find $infile!\n"
		unless -e $infile;
	}


