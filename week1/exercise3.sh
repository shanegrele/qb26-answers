#!/bin/bash

#Configure Table Browser to Assembly hg19, Group Regulation, Table NHEK, Region Genome, Output format BED, Output filename nhek.bed
#Subset regions classified as 1_Active and 12_Repressed in to separate files
#Repeat for NHLF
#Create four .bed files corresponding to 1_Active and 12_Repressed for NHEK and NHLF
    grep 1_Active nhek.bed > nhek-active.bed
    grep 12_Repressed nhek.bed >nhek-repressed.bed
    grep 1_Active nhlf.bed > nhlf-active.bed
    grep 12_Repressed nhlf.bed > nhlf-repressed.bed

#Confirm that your files have 14013, 32314, 14888, and 34469 lines
    wc -l nhek-active.bed #14013 lines
    wc -l nhek-repressed.bed #32314 lines
    wc -l nhlf-active.bed #14888 lines
    wc -l nhlf-repressed.bed #34469 lines

#Construct a bedtools command to test where there is any overlap between 1_Active and 12_Repressed in a given condition (aka mutually exclusive)
    bedtools intersect -a nhek-active.bed -b nhek-repressed.bed | wc -l #0 overlap
    bedtools intersect -a nhlf-active.bed -b nhlf-repressed.bed | wc -l #0 overlap

#Construct two bedtools intersect [OPTIONS] -a nhek-active.bed -b nhlf-active.bed commands, one to find regions that are active in NHEK and NHLF, and one to find regions that are active in NHEK but not active in NHLF
    #active in both:
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed | wc -l #12174 lines
    
    #active in NHEK but not in NHLF: (-v reports entries only in A not B)
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -v | wc -l #2405 lines

#How many features are output by the first command? by the second command?
#Do these two numbers add up to the original number of lines in nhek-active.bed?
    #12174 and 2405, these numbers add to 14579 lines, not the 14013 found

#If not, how can you adjust your first command to only report one feature per overlap?
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -u | wc -l 
    #11608 lines, this adds with 2405 to the correct 14013, -u stops duplicates present in A from being reported as a new line in B

#Construct three bedtools intersect commands to see the effect of using the arguments -f 1, -F 1, and -f 1 -F 1 when comparing -a nhek-active.bed -b nhlf-active.bed
#Visualize the first result for each three commands by pasting the coordinates into UCSC Genome Browser e.g. chr1 25558413 25559413
#Click on Zoom out 3x to get a better perspective.
#How does the relationship between the NHEK and NHLF chromatin state change as you alter the overlap parameter?
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -f 1 | head -1
    #chr1	25558413 25559413, the NHEK active(red) feature is shorter than the NHLF one

    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -F 1 | head -1
    #chr1	19923013 19924213, the NHEK active(red) feature is longer than the NHLF one
    
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -f 1 -F 1 | head -1
    #chr1	1051137	1051537, the NHEK active(red) features and the NHLF ones are of equal length

#Construct three bedtools intersect commands to identify the following types of regions. Use UCSC Genome Browser to save one PDF image for each of the three types of regions. Describe the chromatin state across all nine conditions.
    #Active in NHEK, Active in NHLF
    bedtools intersect -a nhek-active.bed -b nhlf-active.bed -f 1 -F 1 | head -1 
    #chr1	1051137	1051537, this region is active in the other cell lines as well

    #Active in NHEK, Repressed in NHLF
    bedtools intersect -a nhek-active.bed -b nhlf-repressed.bed | head -1 
    #chr1	1981140	1981540, this region is active in some cell lines, repressed in others, and weakly active in others

    #Repressed in NHEK, Repressed in NHLF
    bedtools intersect -a nhek-repressed.bed -b nhlf-repressed.bed -f 1 -F 1 | head -1 
    #chr1	238137	242737, this region is repressed in the other cell lines as well