args<-commandArgs(TRUE);

library("RColorBrewer")

exonHM <- as.matrix(read.table(args[1], header = FALSE))
exonHM[lower.tri(exonHM)] <- NA
pdf (args[2]);
heatmap(log(exonHM+2), Rowv = NA, Colv = NA, labRow=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"), labCol=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"), margins = c(8, 8), revC=TRUE, col=brewer.pal(9,"YlOrRd"), scale="none")
dev.off();

avgExonHM <- as.matrix(read.table(args[3], header = FALSE))
avgExonHM[lower.tri(avgExonHM)] <- NA
pdf (args[4]);
heatmap(log(avgExonHM+1), Rowv = NA, Colv = NA, labRow=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"), labCol=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"), margins = c(8, 8), revC=TRUE, col=brewer.pal(9,"YlGn"), scale="none")
dev.off();

pcHM <- as.matrix(read.table(args[5], header = FALSE))
pcHM[lower.tri(pcHM)] <- NA
pdf (args[6]);
heatmap(log(pcHM+2), Rowv = NA, Colv = NA, labRow=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"), labCol=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"), margins = c(8, 8), revC=TRUE, col=brewer.pal(9,"YlOrRd"), scale="none")
dev.off();

avgPCHM <- as.matrix(read.table(args[7], header = FALSE))
avgPCHM[lower.tri(avgPCHM)] <- NA
pdf (args[8]);
heatmap(log(avgPCHM+1), Rowv = NA, Colv = NA, labRow=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"), labCol=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"), margins = c(8, 8), revC=TRUE, col=brewer.pal(9,"YlGn"), scale="none")
dev.off();
