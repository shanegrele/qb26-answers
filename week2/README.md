#This will be the README for my week2 answers.

#Preparation
##confirm that the dataset is already on the computer 
ls -l BYxRm
##the BYxRM dataset needed to be unzipped:
tar xvf BYxRM.tar 

##skim the genomes 
cd BYxRM
% ls
less BYxRM_GenoData.txt 

#Exercise 1
##1.1

#variants % ls -lh
#the original fq.gz file for A01_09 is 42.7MB, I checked the file size in Data/BYxRM/fastq/A01_09.fq.gz.
47M A01_09.bam
32K A01_09.bam.bai
166M A01_09.sam

The .sam file is larger than the fq.gz file because it is uncompressed. The .sam file also has things like name, quality score, etc.
The .bam file is smaller than the .sam file beacuse it is a compressed binary version of the .sam file.

##1.2: how the for loop would run for A01_23 (the 3rd item)
bwa mem -t 4 -R @RG\tID:A01_23\tSM:A01_23 ../genomes/sacCer3.fa /Users/cmdb/Data/BYxRM/fastq/A01_23.fq.gz
