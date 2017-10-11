#clean gff3 (changed extension from gff3 to tsv)
require(dplyr)
require(tidyr)
colnames(data)<-c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
data <- read.table(file = 'Pvulgaris_442_v2.1.gene.tsv', sep = '\t', header = FALSE, skip = 3)
#extract genes only 
data <- filter(data, type=="gene")
#sort by chromosome/scaffold and by starting coordinate
data <- data[with(data, order(seqid,start)), ]
write.table(test, file='test.tsv', quote=FALSE, sep='\t', col.names = NA)