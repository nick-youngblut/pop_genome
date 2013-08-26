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

my ($verbose);
my ($db_in, $tree_meta_in, $rec_in, $sum_in);
my $ranger_runID;
GetOptions(
	   "database=s" => \$db_in,
	   "metadata=s" => \$tree_meta_in, 
	   "runID=s" => \$ranger_runID,
	   "reconcilliation=s" => \$rec_in,
	   "summary=s" => \$sum_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error

### MAIN
my $sql_r = get_sql();
make_db($sql_r, $ARGV[0], \@tables);


### Subroutines

### Subroutines

	


__END__

=pod

=head1 NAME

ranger_dtl_parse.pl -- parse ranger-dtl output into a *txt table

=head1 SYNOPSIS

ranger_dtl_parse.pl < ranger.dtl.out

=head2 options

=over

=item -prefix

Output file prefix. [ranger-dtl_parse]

=item -h	This help message

=back

=head2 For more information:

perldoc ranger_dtl_parse.pl

=head1 DESCRIPTION

The output are 2 tab-delimited tables: "*_rec.txt" & "*summary.txt"

=head2 Columns of "*_rec.txt" output

=over

=item * 	Tree ID (order in newick input to ranger-dtl)

=item * 	Gene tree node name

=item * 	Species tree node name

=item * 	Category (Duplication, Transfer, Speciation, or Leaf)

=item * 	Recipient (species tree node name)

=back

=head2 Columns of "*_summary.txt" output

=over

=item * 	Tree ID (order in newick input to ranger-dtl)

=item * 	Minimum reconciliation cost

=item * 	Duplication count

=item * 	Transfer count

=item * 	Loss count

=back

=head1 EXAMPLES

=head2 Usage: 

ranger_dtl_parse.pl < ranger-dtl.out 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

