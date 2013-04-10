#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'ITEP_PopGen' ) || print "Bail out!\n";
}

diag( "Testing ITEP_PopGen $ITEP_PopGen::VERSION, Perl $], $^X" );
