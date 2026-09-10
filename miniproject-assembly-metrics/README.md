# Mini Project Assembly Metrics

## Compare C. remanei assemblies

https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA248909/caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa.gz

https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA248911/caenorhabditis_remanei.PRJNA248911.WBPS19.genomic.fa.gz

https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA53967/caenorhabditis_remanei.PRJNA53967.WBPS19.genomic.fa.gz

https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis_remanei/PRJNA577507/caenorhabditis_remanei.PRJNA577507.WBPS19.genomic.fa.gz

## uncompressed file sizes

115M  caenorhabditis_remanei.PRJNA248909.WBPS19.genomic.fa

121M  caenorhabditis_remanei.PRJNA248911.WBPS19.genomic.fa

141M  caenorhabditis_remanei.PRJNA53967.WBPS19.genomic.fa

127M  caenorhabditis_remanei.PRJNA577507.WBPS19.genomic.fa


## python script for contig#, contig length, avg# of contigs
added the python script for this to my github with git add, git commit

to run the script, go to the command line, make sure you're in the correct directory, and type ./assembly-metrics.py followed by the file name

the answers are:

909:
1591,
118549266,
74512.42363293526

911:
912,
124541912,
136559.11403508772

967:
3670,
145442736,
39630.17329700272

507:
187,
130480874,
697758.6844919786

## trend: the older files have more contigs, and a lower average length per contig, likely indicating lower quality