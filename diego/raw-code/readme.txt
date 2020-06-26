This folder contains the three files with the code useful for processing metagenome raw reads. The name "raw" code means that the user
need to provide names/locations of the files that are going to be used. 
The files are:
  quality-trimming.sh- Code that will write two scripts. The first one is to check the quality of the reads; the second one is for trimming the reads.
  metaTOmags.sh - This is a code that will generate four scripts to assemble, bin, assing(taxonomic information) and reassemble(in order to contrive MAGs).
  tax-data-manipulation-phyloseq.R - Code made in R that use phyloseq package to manipulate the information from the taxonomic assigner(kraken).
