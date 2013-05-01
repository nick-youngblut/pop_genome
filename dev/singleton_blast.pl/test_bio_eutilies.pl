#!/usr/bin/env perl

use Bio::DB::EUtilities;
use Data::Dumper;
 
my $id      = 527031;
 
my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'taxonomy',
                                       -id    => $id );

#print Dumper $factory;

my ($name)  = $factory->next_DocSum->get_contents_by_name('ScientificName');
 
print "$name\n";

