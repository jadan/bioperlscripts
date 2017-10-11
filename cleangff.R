#clean gff3 (changed extension from gff3 to tsv)
require(dplyr)
require(tidyr)
data <- read.table(file = 'Pvulgaris_442_v2.1.gene.tsv', sep = '\t', header = FALSE, skip = 3)
colnames(data)<-c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#extract genes only 
data <- filter(data, type=="gene")
#sort by chromosome/scaffold and by starting coordinate
data <- data[with(data, order(seqid,start)), ]
write.table(data, file='genes.gff3', quote=FALSE, sep='\t', col.names = F, row.names=F)