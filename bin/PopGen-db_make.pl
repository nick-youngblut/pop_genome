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
#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
GetOptions(
	   "replace" => \$replace,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error
$ARGV[0] = "PopGen.sqlite" unless $ARGV[0];

### MAIN
my $sql_r = get_sql();
make_db($sql_r, $ARGV[0]);

### Subroutines
sub make_db{
	my ($sql_r, $db_name) = @_;
	
	# checking if tables specified exists, deleted if yes, dying if no #
	if(-e $db_name){
		die " ERROR: $db_name already exists! Use '-r' to replace\n" unless $replace;
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

	$sql{"Cluster_meta"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS Cluster_meta;

CREATE TABLE Cluster_meta(
clusterID	TEXT	NOT NUll,
runID	TEXT	NOT NULL,
core_var	TEXT	NOT NULL,
unique(ClusterID, runID)
ON CONFLICT REPLACE
);

HERE

	$sql{"Fst"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS Fst;

CREATE TABLE Fst(
clusterID	TEXT	NOT NUll,
runID	TEXT	NOT NULL,
pop	TEXT	NOT NULL,
value	REAL	NOT NULL,
ci_low	REAL	NOT NULL,
ci_high	REAL	NOT NULL,
unique(ClusterID, runID, pop)
ON CONFLICT REPLACE
);

HERE

	$sql{"seqID"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS seqID;

CREATE TABLE seqID(
clusterID	TEXT	NOT NUll,
runID	TEXT	NOT NULL,
pop	TEXT	NOT NULL,
value	REAL	NOT NULL,
unique(ClusterID, runID, pop)
ON CONFLICT REPLACE
);

HERE

	$sql{"dN_dS"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS dN_dS;

CREATE TABLE dN_dS(
clusterID	TEXT	NOT NUll,
runID	TEXT	NOT NULL,
pop	TEXT	NOT NULL,

unique(ClusterID, runID, pop)
ON CONFLICT REPLACE
);

HERE
	
	return \%sql;
	}

__END__

=pod

=head1 NAME

rdtl-db_make.pl -- make database tables for Ranger-DTL data

=head1 SYNOPSIS

rdtl-db_make.pl [options] [DATABASE_NAME]

=head2 Options

=over

=item -r 	Replace existing database. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc rdtl-db_make.pl

=head1 DESCRIPTION

Make all tables in the rdtl-db sqlite3 database.

The default database name is "rdtl-db.sqlite"

=head1 EXAMPLES

=head2 Usage: 

rdtl-db_make.pl

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

