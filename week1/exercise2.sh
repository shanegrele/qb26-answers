#!/bin/bash

#Compare hg19 gene annotations with hg16 [fn. 1]

#Prepare hg16 files as in Exercise 1 with the following modifications

#Download hg16.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg16/bigZips/hg16.chrom.sizes 
  ls -l
  less hg16.chrom.sizes
  wc -l hg16.chrom.size

#Do not rename chrM, only exclude secondary contigs e.g.
  grep -v _ hg16.chrom.sizes > hg16-main.chrom.sizes
  
#makewindows command
  bedtools makewindows -g hg16-main.chrom.sizes -w 1000000 > hg16-1mb.bed 

#make the known canonical from the .tsv file downloaded in the genome browser
  cut -f1-3,5 hg16-kc.tsv > hg16-kc.bed

#intersect command
  bedtools intersect -c -b hg16-kc.bed -a hg16-1mb.bed > hg16-kc-count.bed
  wc -l hg16-kc-count.bed #3085 lines, a bit different than hg19


#Visualize both hg19 and hg16 gene distributions on the same line plots
#Create a combined dataframe using bind_rows( hg19=dfA, hg16=dfB, .id="assembly" )
#Place both distributions on the same plot using aes( color=___ )
(see R notebook)

#Calculate how many genes are unique to each assembly

#How many genes are in hg19?
  wc -l hg19-kc.bed 
#80309 lines, but 80308 genes because of the header line

#How many genes are in hg19 but not in hg16?
#Use intersect with a 1-letter option to find genes with no overlaps
  bedtools intersect -a hg19-kc.bed -b hg16-kc.bed -v > unique_hg19-kc.bed
  wc -l unique_hg19-kc.bed 
#there are 43217 genes present in hg19 but not hg16

#Why are some genes in hg19 but not in hg16?
#hg19 represents an improvement to sequencing technologies, and there are likely more transcripts reflected.

#Answer the same three questions but with respect to hg16

  wc -l hg16-kc.bed
#21365 lines, but 21364 genes because of the header line

  bedtools intersect -a hg16-kc.bed -b hg19-kc.bed -v > unique_hg16-kc.bed
  wc -l unique_hg16-kc.bed
#there are 3458 genes present in hg16 but not hg19

#there might be lines present in hg16 but not hg19 beause the older assemblies were probably mismapped due to issues with mapping repetitive elements.
