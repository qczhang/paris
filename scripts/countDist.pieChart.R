args<-commandArgs(TRUE);

count <- read.table ( args[1], header=TRUE, sep="\t");
pdf ( args[2], height=8, width=12);

lbls = paste (count$type, ":", count$count, sep="" );
pie ( count$count, labels = lbls, main = "Distribution of stem duplex in different categories of RNAs" );

dev.off();

