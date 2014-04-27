# pop_genome

## DESCRIPTION

A collection of scripts used for various population genomics analyses of bacterial/archaeal genomes.
Many build upon ITEP, which is a toolset for pangenome analysis. 

## DEPENDENCIES

Depending on the script:

* BioPerl 
* ITEP: https://github.com/mattb112885/clusterDbAnalysis
* SNAP: http://www.hiv.lanl.gov/content/sequence/SNAP/SNAP.html
* PAML: http://abacus.gene.ucl.ac.uk/software/paml.html
* Arlecore: http://cmpg.unibe.ch/software/arlequin35/
* Mowgli: http://www.atgc-montpellier.fr/Mowgli/
* Ranger-DTL: http://compbio.mit.edu/ranger-dtl/
* Mothur: http://www.mothur.org/wiki/Download_mothur
* PhiPack: http://www.maths.otago.ac.nz/~dbryant/software.html
* sqlite


## INSTALLATION

To install this module, run the following commands:
	
	perl Build.PL
	./Build
	./Build test
	./Build install

## Scripts used in 'Genomic and phenotypic differentiation among Methanosarcina mazei populations from Columbia River sediment'

This repo also includes scripts for other analyses no included in the manuscript.
All the scripts should have documentation (use perldoc) and all work, at least
with my dataset.

Scripts named as "ITEP_*" require ITEP (or files produced by ITEP) to run.

### Calculating Fst, mean sequence identity, & dN/dS for each ITEP sequence cluster

arlecore_Fst_batch.pl

seqID_byPop.pl

SNAP_batch.pl

### Batch runs of Mowgli

Mowgli_batch.pl

### Parsing output from the Quartet Decomposition Server (http://csbl1.bmb.uga.edu/QD/)

quartetDecomp_byPop.pl

quartetDecomp_distCut.pl

### Determine how far a gene is to the end of a contig

ITEP_geneInfoDist2ContigEnd.pl

### Get distance between genes features

ITEP_geneInfo2scaffoldDist.pl

### Simple pipeline for aligning (mafft) and inferring the phylogenies (RAxML) of ITEP gene clusters.

ITEP_align-tree.pl

### Convert an ITEP geneInfo table to a gff

ITEP_geneInfo2gff.pl

### Make a database of output from multiple Ranger-DTL runs

ranger-dtl_parse.pl

ranger-dtl_parse_ITOL.pl

rdtl-db_*

### Calculate feature density for a Circos feature file

circos_featDensity.pl


## SUPPORT AND DOCUMENTATION

After installing, you can find documentation for this module with the
perldoc command.

    perldoc ITEP_PopGen

You can also look for information at:

    RT, CPAN's request tracker (report bugs here)
        http://rt.cpan.org/NoAuth/Bugs.html?Dist=ITEP_PopGen

    AnnoCPAN, Annotated CPAN documentation
        http://annocpan.org/dist/ITEP_PopGen

    CPAN Ratings
        http://cpanratings.perl.org/d/ITEP_PopGen

    Search CPAN
        http://search.cpan.org/dist/ITEP_PopGen/


## LICENSE AND COPYRIGHT

Copyright (C) 2013 Nick Youngblut

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.

