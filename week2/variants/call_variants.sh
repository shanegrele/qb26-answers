#!/bin/bash

# Step 3.1: make the .txt file with the list of .bam files in 'variants'
ls *.bam > bamListFile.txt
# run FreeBayes to discover variants
#first, I ran 
mamba install -y -c conda-forge -c bioconda freebayes vcflib

#then, create unfiltered VCF:
freebayes -f "../genomes/sacCer3.fa" -L bamListFile.txt --genotype-qualities -p 1 > unfiltered.vcf

# the resulting VCF file is unfiltered, meaning that it contains low-confidence calls and also has
# some quirky formatting, so the following steps use a software suite called vcflib to clean it up
# filter the variants based on their quality score and remove sites where any sample had missing data
vcffilter -f "QUAL > 20" -f "AN > 9" unfiltered.vcf > filtered.vcf

# FreeBayes has a quirk where it sometimes records haplotypes rather than individual variants;
# we want to override this behavior
vcfallelicprimitives -kg filtered.vcf > decomposed.vcf

# in very rare cases, a single site may have more than two alleles detected in your sample; while
# these cases may be interesting, they may also reflect technical errors and pose a challenge for
# parsing the data, so we remove them
vcfbreakmulti decomposed.vcf > biallelic.vcf

#this was taking a very long time in the terminal, so I used the downloaded biallelic.vcf version from the bxlab repository!


