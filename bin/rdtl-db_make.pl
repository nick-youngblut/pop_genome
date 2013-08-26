#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
GetOptions(
	   "replace" => \$replace,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error

### MAIN
$sql_r = get_sql();
make_db($sql_r, $ARGV[0], \@tables);

### Subroutines
sub make_db{
	my ($sql_r, $db_name) = @_;
	
	# checking if tables specified exists, deleted if yes, dying if no #
	if(-e $db_name){
		die " ERROR: $db_name already exists! Use '-r' to replace\n" unless $replace
		unlink $db_name;
		}
	
	# adding tables #
	foreach my $table (keys %$sql_r){
		open PIPE, "| sqlite3 $db_name" or die $!;
		print PIPE "BEGIN TRANSACTION;\n";
		print PIPE $sql_r->{$table}; 				# making table
		print PIPE "COMMIT;\n";
		close PIPE;
		}
		
	print STDERR "...sqlite3 database tables created\n";
	}

sub get_sql{
	my %sql; 		# all tables individually 

	$sql{"loci"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS Loci;

CREATE TABLE Loci (
Locus_ID	INTEGER	PRIMARY KEY,
Taxon_ID	TEXT	NOT NULL,
Taxon_Name	TEXT	NOT NULL,
Subtype	TEXT,
Scaffold	TEXT	NOT NULL,
Locus_Start	INTEGER	NOT NULL,
Locus_End	INTEGER	NOT NULL,
Operon_Start	INTEGER,
Operon_End	INTEGER,
CRISPR_Array_Start	INTEGER,
CRISPR_Array_End	INTEGER,
Operon_Status	TEXT	NOT NULL,
CRISPR_Array_Status	TEXT	NOT NULL,
Genbank	TEXT	NOT NULL,
Array_File	TEXT,
Scaffold_count	INTEGER,
File_Creation_Date	TEXT,
Author	TEXT	NOT NULL,
UNIQUE (Taxon_ID, Taxon_name, Scaffold, Locus_Start, Locus_End)
ON CONFLICT IGNORE
);

HERE
	
	return \%sql;
	}

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

