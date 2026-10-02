#!/usr/bin/env python3

af = open("AF.txt", "w")
gt_long = open("gt_long.txt", "w")

for line in open("/Users/cmdb/qb26-answers/week2/variants/biallelic.vcf"):
    #Skip the header lines
    if line.startswith('#'):
        continue
    #Remove the newline character and tab-separate
    fields = line.rstrip('\n').split('\t')
    #Ignore chrM
    if fields[0] =="chrM": continue
    #label the fields
    chrom = fields[0]
    pos = fields[1]
    id = fields[2]
    ref = fields[3]
    alt= fields[4]
    qual= fields[5]
    filter= fields[6]
    info=fields[7]
    format= fields[8]
    sample_fields = fields[9:]
    sample_names = ["A01_62", "A01_39", "A01_63", "A01_35", "A01_31",
                 "A01_27", "A01_24", "A01_23", "A01_11", "A01_09"]


    #get AF from the INFO field, split to make readable
    allele_freq = info.split("AF=")[1].split(";")[0]

    #Write a file into 'af' to get the allele frequency for each variant
    af.write(chrom + "\t" + pos + "\t" + allele_freq + "\n")

    #Write a file into 'gt_long' to get the genotype of each variant
    for sample_name, sample_field in zip(sample_names, sample_fields):
        genotype = sample_field.split(":")[0]
        gt_long.write(chrom + "\t" + pos + "\t" + sample_name + "\t" + genotype + "\n") 
   
#close
af.close()
gt_long.close()