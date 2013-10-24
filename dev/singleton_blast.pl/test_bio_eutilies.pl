#!/usr/bin/env perl

use Bio::DB::EUtilities;
use Bio::DB::Taxonomy;
use Data::Dumper;
 
my $id      = 527031;

# 823;563193;999416

my $db = new Bio::DB::Taxonomy(-source => 'entrez');
my $node= $db->get_Taxonomy_Node(-taxonid => $id);

# get basal node #
my $anc = $node;
while(1){
	last unless $anc->ancestor();
	$anc = $anc->ancestor();
	print "id is ", $anc->id, "\n"; # NCBI taxa id
    print "rank is ", $anc->rank, "\n"; # e.g. species
	print "scientific name is ", $anc->scientific_name, "\n"; # scientific name
	}

#print $node->id, " ", $node->scientific_name, " ", $node->rank, "\n";

my @desc = $db->get_all_Descendents($node);

for my $child ( @extant_children ) {
    print "id is ", $child->id, "\n"; # NCBI taxa id
    print "rank is ", $child->rank, "\n"; # e.g. species
    print "scientific name is ", $child->scientific_name, "\n"; # scientific name
}

exit;


## Bio::DB::EUtilities ##
my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                       -email => 'mymail@foo.bar',
                                       -db    => 'taxonomy',
                                       -id    => $id );

#print Dumper $factory; exit;

#print Dumper $factory->next_DocSum->get_contents_by_name('Status');

for my $doc ($factory->next_DocSum()){
	#print Dumper $doc; exit;
	for my $name ($doc->get_all_names()){
		print join("\t", $name, $doc->get_contents_by_name($name)), "\n";
		}
	}
 
 
#print "$name\n";

