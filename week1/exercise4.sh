#!/bin/bash

#Obtain snps-chr1.bed with only the chromosome 1 Common SNPs
#Configure Table Browser to Assembly hg19, Group Variation, Track Common SNPs(151), Region Position chr1, Output format BED, Output filename snps-chr1.bed
#Confirm that your file has 1091148 lines
    wc -l snps-chr1.bed #1091148 lines

#Use bedtools intersect and hg19-kc.bed to determine which gene has the most SNPs
    bedtools intersect -a hg19-kc.bed -b snps-chr1.bed | sort -k5,5nr | head -1
    #chr1	245912648	246670581	ENST00000490107.6_7	5445

#Visualize the gene using UCSC Genome Browser and click on the gene to see more details
#Describe the gene: 
    #What is the systematic name: ENST00000490107.6_7
    #human readable name: SMYD3
    #position: hg19 chr1:245,912,648-246,670,581
    #size: 668,068
    #exon count: 12
    #Why do you think this gene has the most SNPs?
        #SMYD3 is a large gene, that therefore has more opportunities for mutations.
        #It also has repetitive sequences with high GC-content, leading to variation in copy number.

#Determine which SNPs lie within vs outside of a gene

    #Create a subset of SNPs using bedtools sample -n 20 -seed 42 (20 random samples)
    bedtools sample -n 20 -seed 42 -i snps-chr1.bed > subset.bed | head -5
    #chr1	247118389	247118390	rs61852577	0	+
    #chr1	11839540	11839541	rs200540772	0	+
    #chr1	36733538	36733539	rs12726228	0	+
    #chr1	201437831	201437832	rs35383942	0	+
    #chr1	174840899	174840900	rs1883139	0	-

    #Use bedtools sort to sort the subset of SNPs
    bedtools sort -i subset.bed > subset-sorted.bed | head -5
    #chr1	3810505	    3810506	    rs78397137	0	+
    #chr1	11638083	11638084	rs6698664	0	+
    #chr1	11839540	11839541	rs200540772	0	+
    #chr1	19020850	19020851	rs71645417	0	+
    #chr1	19821839	19821840	rs2088825	0	+

    #Use bedtools sort to sort hg19-kc.bed
    bedtools sort -i hg19-kc.bed > sorted-hg19.bed | head -5
    #chr1	10369	10582	ENST00000833856.1_2
    #chr1	11425	14409	ENST00000832828.1_1
    #chr1	12009	13670	ENST00000450305.2_3
    #chr1	14360	29367	ENST00000831158.1_1
    #chr1	14695	24886	ENST00000488147.2_6

    #Use bedtools closest -d on the two sorted files, with -t first to break ties
    bedtools closest -d -t first -a subset-sorted.bed -b sorted-hg19.bed > out.bed | head -5
    #chr1	3810505	    3810506	    rs78397137	0	+	chr1	3805696	    3816836	    ENST00000361605.4_7	0
    #chr1	11638083	11638084	rs6698664	0	+	chr1	11653741	11655507	ENST00000793460.1_2	15658
    #chr1	11839540	11839541	rs200540772	0	+	chr1	11822249	11849642	ENST00000688073.1_8	0
    #chr1	19020850	19020851	rs71645417	0	+	chr1	18957339	19075360	ENST00000420770.7_5	0
    #chr1	19821839	19821840	rs2088825	0	+	chr1	19823503	19890741	ENST00000816783.1_2	1664

    #How many SNPs are inside of a gene?
    awk '$11 == 0' out.bed | wc -l #pull any that have a column 11(distance)as zero, and count the # of lines
    #5

    #What is the range of distances for the ones outside a gene?
    cut -f 11 out.bed | sort -n
    #0
    #0
    #0
    #0
    #0
    #1664 <- shortest distance: 1664bp
    #4407
    #15658
    #11716456
    #16415963
    #23280856
    #25993140
    #44313390
    #116776105
    #125650729
    #135515456
    #135625913
    #150277473
    #162112388
    #207792946 <- longest distance: 207792946bp