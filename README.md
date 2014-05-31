# pop_genome

## DESCRIPTION

A collection of scripts used for comparative population genomics
analyses of large microbial genome datasets.
Many of the scripts are build upon[ITEP](https://github.com/mattb112885/clusterDbAnalysis), 
which is a toolset for pangenome analysis. 

## DEPENDENCIES

Depending on the script:

* [BioPerl](http://www.bioperl.org/wiki/Main_Page)
* [ITEP](https://github.com/mattb112885/clusterDbAnalysis)
* [SNAP](http://www.hiv.lanl.gov/content/sequence/SNAP/SNAP.html)
* [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
* [Arlecore](http://cmpg.unibe.ch/software/arlequin35/)
* [Mowgli](http://www.atgc-montpellier.fr/Mowgli/)
* [Ranger-DTL](http://compbio.mit.edu/ranger-dtl/)
* [Mothur](http://www.mothur.org/wiki/Download_mothur)
* [PhiPack](http://www.maths.otago.ac.nz/~dbryant/software.html)
* sqlite


## INSTALLATION

   git clone https://github.com/nyoungb2/pop_genome.git

To install this module, run the following commands:
	
	perl Build.PL
	./Build
	./Build test
	./Build install

__Or__ just add the './bin/' folder to your $PATH.


## Getting help

### Script documentation

    SCRIPT.pl -h
    perldoc SCRIPT.pl

### Post issues on GitHub

    [Issues](https://github.com/nyoungb2/pop_genome/issues)


## ITEP Extensions

Scripts named as "ITEP_*" require ITEP (or files produced by ITEP) to run.

### Some features

* Convert the ITEP geneInfo table to a gff3 file:
 * ITEP_geneInfo2gff.pl
* Get info on how far apart genes are on each genome
 * ITEP_geneInfo2scaffoldDist.pl
* For draft genes: get the distance of a gene to a contig end 
(helps for finding artificial gene trucations)
 * ITEP_geneInfoDist2ContigEnd.pl
* Make a ITOL compatible metadata table of gene copies per taxon-gene_cluster
 * ITEP_geneInfo2ITOL.pl
* My own alignment and phylogeny inference pipeline to use with ITEP
 * ITEP_align-tree.pl
* Add COG info to a table containing PED IDs (eg., the geneInfo table)
 * ITEP_clusterAddCOG.pl


## Other Scripts 

These scripts provide convenient methods for using many population genomics
software on large datasets. 

Many of these were used for the analysis in:
'N.D. Youngblut, J.S. Wirth, J.R. Henriksen , W.W. Metcalf, R.J. Whitaker
_Genomic and phenotypic differentiation among Methanosarcina mazei populations from Columbia River sediment_,
_in prep_'



### Calculating Fst, mean sequence identity, & dN/dS for each ITEP sequence cluster

* arlecore_Fst_batch.pl
* seqID_byPop.pl
* SNAP_batch.pl

### Batch runs of Mowgli

* Mowgli_batch.pl

### Parsing output from the Quartet Decomposition Server (http://csbl1.bmb.uga.edu/QD/)

* quartetDecomp_byPop.pl
* quartetDecomp_distCut.pl


### Make a database of output from multiple Ranger-DTL runs

* ranger-dtl_parse.pl
* ranger-dtl_parse_ITOL.pl
* rdtl-db_*

### Calculate feature density for a Circos feature file

* circos_featDensity.pl



## LICENSE AND COPYRIGHT

Copyright (C) 2013 Nick Youngblut

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.

