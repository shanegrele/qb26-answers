# This will be the README for my week2 answers.

# Preparation
## Confirm that the dataset is already on the computer 
ls -l BYxRm
## The BYxRM dataset needed to be unzipped:
tar xvf BYxRM.tar 
## Skim the genomes 
cd BYxRM

% ls

less BYxRM_GenoData.txt 

# Exercise 1
## Commands 

**see map_reads.sh**

## Question 1.1

the original fq.gz file for A01_09 is 42.7MB, I checked the file size in Data/BYxRM/fastq/A01_09.fq.gz.

47M A01_09.bam

32K A01_09.bam.bai

166M A01_09.sam

The .sam file is larger than the .fq.gz file because it is uncompressed. The .sam file also has things like name, quality score, etc.
The .bam file is smaller than the .sam file beacuse it is a compressed binary version of the .sam file.

## Question 1.2: how the for loop would run for A01_23 (the 3rd item)
bwa mem -t 4 -R @RG\tID:A01_23\tSM:A01_23 ../genomes/sacCer3.fa /Users/cmdb/Data/BYxRM/fastq/A01_23.fq.gz

# Exercise 2

## Commands
**step 2.1:** 

samtools view -H A01_09.bam

samtools view A01_09.bam | head -n 3

  output:
  
  HWI-ST387_0114:5:47:16737:7355#0	0	chrI	29	60	76M	*	0	0	ACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCA	@FGDFGFGGDGGGGHHGHEHHHGHHHGFHHFFDC>DEEEEDDB< BE>?7?>C@@@@CDBDDF?E2FBDAA>361::	NM:i:0	MD:Z:76	AS:i:76	XS:i:22	RG:Z:A01_09
  
  HWI-ST387_0114:5:23:10710:84112#0	0	chrI	47	60	76M	*	0	0	CACCCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTAC	BE8E=EGGFGFEEFEFGGGFCFFFFD?BDDDDFEBD:A?A:ED=BDA?ADD28..8:>6>;6>??BFDED.DAEBE	NM:i:0	MD:Z:76	AS:i:76	XS:i:0	RG:Z:A01_09
  
  HWI-ST387_0114:5:42:14962:12601#0	16	chrI	47	60	4S72M	*	0	0	ACAGCACCCACACACACACATCCTAACACTACCCTAACACAGTCCTAATCTAACCCTGGCCAACCTGTCTCTCAAC	###############BB>B@.CAAC<9;@61=A>/BADAD:E;FECDFCB=FCFGFAGFGG@DFGEFFFDEDGGGE	NM:i:1	MD:Z:38C33	AS:i:67	XS:i:26	RG:Z:A01_09

**step 2.2**

samtools flagstat A01_09.bam > A01_09.flagstat

cat A01_09.flagstat

**step 2.3**

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 10

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 12

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 24

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 25

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 28

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 32

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 36

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 40

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 63

grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 64


## Question 2.1 
There are 17 @SQ lines, and they are each one sequence in the reference genome that the reads were aligned to. 
It gives you the name of the sequence(SN) and the length(LN) in bp.

## Question 2.2
The 1st alignment is aligned to chromosome 1 at position 29. The CIGAR score, 76M, means that there are 76 bases that match. Other CIGAR strings can indicate insertions, deletions, etc.

## Question 2.3 
@RG and RG:Z comes from the ID and SM that was attributed to the samples when indexing.
The ID in @RG is used as the batch for sequencing.
The SM in RG:Z is the sample name. This specifies which individual or organism the reads come from.
The format of SAM requires both, and even though our files are 1:1, this is useful to keep track of samples.

## Question 2.4
100% of the reads mapped, which could be feasible if the genome isn't very big or there isn't great coverage, but also could be suspicious because we'd expect at least a small amount of un-mapped reads that are lower quality, or maybe some differences across yeast samples.

## Question 2.5
The zeroes are saying that criteria for paired-end reads does not apply. Each fragment was sequenced with one end only, and reads do not have a 'mate'. Sequencing was single-end and not paired-end.

## Question 2.6
A01_09, AO1_24, A01_31, A01_34, A01_62, and A01_63 are mostly grey, and likely carry BY ancestry at that region.
A01_11, A01_23, A01_27, A01_35, A01_39, are more colorful, and likely carry RM ancestry at that region.

I then ran this in the terminal for all the samples to look at the chromosome regions,

example for A01_09: 

input: grep -F -e 'chr01_27915' -e 'chr01_28323' -e 'chr01_28652' -e 'chr01_29667' ~/Data/BYxRM/BYxRM_GenoData.txt | cut -f 10

output: 
B
B
B
B

I found that all but A01_39 agree with my visual call.

# Exercise 3

## Commands 

**see call_variants.sh**

## Question 3.1

**ran less -S biallelic.vcf**

The last 10 columns in the #CHROM line come from the sample names, which are from the name of the .bam files.


## Question 3.2

With the command that was run to get the vcf file originally with freebayes, the "-p 1" is 1-ploid, or haploid.
If it were diploid, or "-p 2", then the resultant genotypes could be heterozygous, and freebayes would have to call whether the sample matches the REF or the ALT as 0/0, 0/1, or 1/1, instead of the haploid 0 or 1.













