args<-commandArgs(TRUE);

sizeDist <- read.table (args[1], header=TRUE);
attach(sizeDist);

pdf (args[2]);
library ( ggplot2 );
ggplot(sizeDist, aes(type, log10(genomeSize))) + geom_violin(colour="black",fill="red") + scale_y_continuous( breaks = seq(log10(min(sizeDist$genomeSize)),log10(max(sizeDist$genomeSize)),len=5), labels = round(10^(seq(log10(min(sizeDist$genomeSize)),log10(max(sizeDist$genomeSize)),len=5))))

dev.off();

detach(sizeDist);
