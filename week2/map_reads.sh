!#/bin/bash

#Step 1.1: Prepare the reference sample
    #first, create a directory called genomes in week2
    week2 % mkdir genomes
    #then copy and unzip the yeast genome files into the genomes directory
    cd genomes
    cp ~/Data/References/sacCer3/sacCer3.fa.gz .
    gunzip sacCer3.fa.gz
    bwa index sacCer3.fa

    #check with 
    ls -l #sacCer3.fa	sacCer3.fa.amb	sacCer3.fa.ann	sacCer3.fa.bwt	sacCer3.fa.pac	sacCer3.fa.sa 


#Step 1.2: Map one sample
    #first, create a directory called variants in week2
    week2 % mkdir variants

    #then, you can run this command
    cd .. /variants
    bwa mem -t 4 -R "@RG\tID:A01_09\tSM:A01_09" ../genomes/sacCer3.fa ~/Data/BYxRM/fastq/A01_09.fq.gz > A01_09.sam

#Step 1.3
    samtools sort -@ 4 -O bam -o A01_09.bam A01_09.sam
    samtools index A01_09.bam

#compare file sizes (see README)
    variants % ls -lh
    #the original fq.gz file for A01_09 is 42.7MB

#Step 1.4
#practice echo
    my_sample=A01_09
    echo ${my_sample} #this prints A1_09

#practice echo
    for my_sample in A01_09 A01_11 A01_23
    for> do
    for>    echo "Now processing" ${my_sample}
    for> done
    #Now processing A01_09
    #Now processing A01_11
    #Now processing A01_23

#Step 1.5: the 10 segregants we will analyze
    for my_sample in A01_09 A01_11 A01_23 A01_24 A01_27 A01_31 A01_35 A01_39 A01_62 A01_63
    do
        echo "***" ${my_sample}

        # align reads to the reference genome
        bwa mem -t 4 -R "@RG\tID:${my_sample}\tSM:${my_sample}" ../genomes/sacCer3.fa ~/Data/BYxRM/fastq/${my_sample}.fq.gz > ${my_sample}.sam

        # sort the alignments by position and convert to BAM
        samtools sort -@ 4 -O bam -o ${my_sample}.bam ${my_sample}.sam

        # index the BAM file
        samtools index ${my_sample}.bam
    done
#(see README for question 1.2)