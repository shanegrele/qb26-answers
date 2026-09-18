#!/bin/bash

#makewindows command
bedtools makewindows -g hg19-main.chrom.sizes -w 1000000 > hg19-1mb.bed 

#intersect command
bedtools intersect -c -b hg19-kc.bed -a hg19-1mb.bed > hg19-kc-count.bed


