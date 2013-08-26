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
$ARGV[0] = "rdtl-db.sqlite" unless $ARGV[0];

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

	$sql{"Tree_meta"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS Tree_meta;

CREATE TABLE Tree_meta(
TreeID	TEXT	NOT NULL,
RunID	TEXT	NOT NULL,
ClusterID	TEXT	NOT NULL,
Date	Date,
UNIQUE (TreeID, RunID)
ON CONFLICT REPLACE
);

HERE

	$sql{"Ranger_run"} = <<HERE;
DROP TABLE IF EXISTS Ranger_run;

CREATE TABLE Ranger_run(
Ranger_runID	TEXT	PRIMARY KEY,
TreeID	TEXT	NOT NULL,
D_cost	INTEGER	NOT NULL,
T_cost	INTEGER NOT NULL,
L_cost 	INTEGER	NOT NULL,
Run_date	Date,
UNIQUE (Ranger_runID, TreeID)
ON CONFLICT REPLACE
);

HERE

	$sql{"Node"} = <<HERE;
DROP TABLE IF EXISTS Node;

CREATE TABLE Node(
Ranger_runID	TEXT	NOT NULL,
TreeID	TEXT	NOT NULL,
Category	TEXT	NOT NULL,
Gene_NodeID	TEXT	NOT NULL,
Species_NodeID	TEXT	NOT NULL,
T_recipient	TEXT,
UNIQUE (Ranger_runID, TreeID, Gene_NodeID)
ON CONFLICT REPLACE
);

HERE

	$sql{"Tree"} = <<HERE;
DROP TABLE IF EXISTS Tree;

CREATE TABLE Tree(
Ranger_runID	TEXT	NOT NULL,
TreeID	TEXT	NOT NULL,
Min_rec_cost	INTEGER	NOT NULL,
D_cnt	INTEGER	NOT NULL,
T_cnt	INTEGER	NOT NULL,
L_cnt	INTEGER	NOT NULL,
UNIQUE (Ranger_runID, TreeID)
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

