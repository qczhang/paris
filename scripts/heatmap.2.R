args<-commandArgs(TRUE);

if (!require("gplots")) { install.packages("gplots", dependencies = TRUE); library(gplots) }
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE); library(RColorBrewer) }

exonHM <- as.matrix(read.table(args[1], header = FALSE))
exonHM[lower.tri(exonHM)] <- NA
avgExonHM <- as.matrix(read.table(args[3], header = FALSE))
avgExonHM[lower.tri(avgExonHM)] <- NA
pcHM <- as.matrix(read.table(args[5], header = FALSE))
pcHM[lower.tri(pcHM)] <- NA
avgPCHM <- as.matrix(read.table(args[7], header = FALSE))
avgPCHM[lower.tri(avgPCHM)] <- NA

#lmat <- rbind(c(0,3),c(2,1),c(0,4))
lhei <- c(1,5)
lwid <- c(1,4)

pdf (args[2]);
par(cex.main=0.90)
heatmap.2(log2(exonHM+2),
        cellnote = exonHM,  # same data set for cell labels
	lhei = lhei,
	lwid = lwid,
        main = "Counts of connections\nwithin/across different exons", # heat map title
        notecol="black",      # change font color of cell labels to black
        density.info="none",  # turns off density plot inside color legend
        trace="none",         # turns off trace lines inside the heat map
        margins =c(8,8),     # widens margins around plot
        col=brewer.pal(9,"YlOrRd"),       # use on color palette defined earlier 
        labRow=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"),
        labCol=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"),
        scale="none",
        Rowv=NA,            # turn off column clustering
        Colv=NA)            # turn off column clustering
dev.off();

pdf (args[4]);
par(cex.main=0.90)
heatmap.2(avgExonHM,
        cellnote = round(avgExonHM, digits = 2),  # same data set for cell labels
	lhei = lhei,
	lwid = lwid,
        main = "Counts of connections within/across\ndifferent exons normalized by size", # heat map title
        notecol="black",      # change font color of cell labels to black
        density.info="none",  # turns off density plot inside color legend
        trace="none",         # turns off trace lines inside the heat map
        margins =c(8,8),     # widens margins around plot
        col=brewer.pal(9,"YlGn"),       # use on color palette defined earlier 
        labRow=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"),
        labCol=c("exon 1", "exon 2", "exon 3", "exon -3", "exon -2", "exon -1"),
        scale="none",
        Rowv=NA,            # turn off column clustering
        Colv=NA)            # turn off column clustering
dev.off();

pdf (args[6]);
par(cex.main=0.90)
heatmap.2(log2(pcHM+2),
        cellnote = pcHM,  # same data set for cell labels
	lhei = lhei,
	lwid = lwid,
        main = "         Counts of connections within/across different regions\n of protein coding genes", # heat map title
        notecol="black",      # change font color of cell labels to black
        density.info="none",  # turns off density plot inside color legend
        trace="none",         # turns off trace lines inside the heat map
        margins =c(8,8),     # widens margins around plot
        col=brewer.pal(9,"YlOrRd"),       # use on color palette defined earlier 
        labRow=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"),
        labCol=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"),
        scale="none",
        Rowv=NA,            # turn off column clustering
        Colv=NA)            # turn off column clustering
dev.off();

pdf (args[8]);
par(cex.main=0.90)
heatmap.2(avgPCHM,
        cellnote = round(avgPCHM, digits = 2),  # same data set for cell labels
	lhei = lhei,
	lwid = lwid,
        main = "          Counts of connections within/across different regions\n        of protein coding genes normalized by size", # heat map title
        notecol="black",      # change font color of cell labels to black
        density.info="none",  # turns off density plot inside color legend
        trace="none",         # turns off trace lines inside the heat map
        margins =c(8,8),     # widens margins around plot
        col=brewer.pal(9,"YlGn"),       # use on color palette defined earlier 
        labRow=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"),
        labCol=c("5UTR", "start codon", "CDS", "stop codon", "3UTR"),
        scale="none",
        Rowv=NA,            # turn off column clustering
        Colv=NA)            # turn off column clustering
dev.off();
