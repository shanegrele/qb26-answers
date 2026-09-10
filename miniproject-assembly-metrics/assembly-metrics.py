#!/usr/bin/env python3

import sys
import fasta

file=open(sys.argv[1])
genomes=fasta.FASTAReader(file)

#counts the number of contigs when you type the file name at the command line
#also adds all the sequence lengths together in the for loop
count=0
total_length=0
for ident, sequence in genomes:
    count += 1
    total_length += int(len(sequence))
print(count)
print(total_length)

#gets the average length as the total length divided by the # of contigs
average_length=total_length/count
print(average_length)

file.close() 

#Count the number of contigs
#Determine each sequence length using len() and sum up the total length
#Print the “Number of contigs: “, “Total length: “, and “Average length: “

#1591
#118549266
#74512.42363293526

#912
#124541912
#136559.11403508772

#3670
#145442736
#39630.17329700272

#187
#130480874
#697758.6844919786