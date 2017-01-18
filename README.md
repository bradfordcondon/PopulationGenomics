
#Overview

1-17-17
Bradford Condon <bradford.condon@gmail.com>

The collection will include all scripts required to 

* Automate working with ITS sequences
* Align orthologs
* Calculate population genomics statistics (ie. FST)

###Requirements and set up

* [ITSx](http://microbiology.se/software/itsx/)
* BioPerl (for BLAST-based ITS)



#ITS

ITS regions are important barcodes for fungi.  This repo contains two sets of scripts for extracting them: ITSx and BLAST.  The main problem with an ITS approach is the sequences are often left out of assemblies for next-generation sequenced genomes.  Because of this, tools that extract ITS regions are less helpful than an approach which amplifies and sequences ITS regions from the beginning. 


###ITSx based

* **run_ITS.pl**

Runs ITSx on all FASTA files in the current directory, or, specificy a directory with -w.

* **process_ITSX.pl**

Processes ITSx output.  Peforms filtering (regions > 75, < 400), outputs two FASTA files- one for ITS1, one for ITS2.


###Blast-based

* run_BLAST

Runs BLAST given a FASTA query file an directory of blast databases.	

* combine_ITS_blasts.pl

Combines BLAST output


###Output

After generating the ITS sequences, they can be aligned and trimmed as follows:

```
muscle3 -in data/ITS1_all_Loci_100min.fasta -out data/ITS1_100min_aligned.fasta
Gblocks data/ITS1_100min_aligned.fasta -t=d -b5=n
```

#Population genomics

Given a set of aligned orthologs (each gene in a seperate file)


* cleanNames_orthologs.sh

Simple shell script that renames FASTA headers.  Necessary for how input was named.
```
GRF52 MGG_17959T0_GRF52_scaffold00200.g7712.t1 GRF52_scaffold00200.g7712.t1 
```
Contains gene model & scaffold location.  

###BLAST

In addition to starting with a set of orthologs, I wanted to generate my own set for a specific gene.  With that in mind, I use the two beow scripts to BLAST a query gene against all genomes, and then concatenate the results.

* run_BLAST_for_orthologs.pl
* Read_and_stitch_orthologs.pl
