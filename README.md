
#Overview

1-17-17
Bradford Condon <bradford.condon@gmail.com>

This repo is an early draft.

The collection will include all scripts required to 

* Automate working with ITS sequences
* Align orthologs
* Calculate population genomics statistics (ie. FST)


###Requirements and set up

[ITSx]: http://microbiology.se/software/itsx/

BioPerl (for BLAST-based ITS)

#Scripts


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