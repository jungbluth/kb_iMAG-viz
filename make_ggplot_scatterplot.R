#!/usr/bin/env Rscript

library(ggplot2)
dat1<-read.csv("_Plot-input-data.tsv", sep="\t", header=T)
head(dat1)
p1 <- ggplot(dat1) + theme_bw() + geom_point(aes(x=dat1$PC1, y=dat1$PC2, alpha=1.0, fill=Order, shape=GenomeType)) + scale_shape_manual(values=c(21, 22, 24, 25)) + guides(alpha=FALSE, fill=TRUE, size=FALSE, color=FALSE, reverse=TRUE) + guides(fill=guide_legend(title="Order")) + xlab("PC1") + ylab("PC2") + labs(shape = "Genome Type", fill = "Order")
#p1 <- p1 + geom_text(aes(x=dat1$PC1,y=dat1$PC2,label=dat1$genomeID,hjust=0, vjust=0)) + theme(text=element_text(color="black",size=9),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=9)) # add labels
#p1 <- ggplot(dat1) + theme_bw() + geom_point(aes(x=dat1$PC1, y=dat1$PC2, alpha=0.5, size=Completeness, fill=Order, shape=GenomeType))  + geom_point(aes(x=dat1$PC1, y=dat1$PC2, alpha=1, stroke=1.0, color=I("black"),fill=Order,shape=GenomeType)) + scale_shape_manual(values=c(21, 22, 24, 25)) + guides(alpha=FALSE, fill=FALSE, size=FALSE, color=FALSE) + xlab("PC1") + ylab("PC2") + labs(shape = "Genome Type", fill = "Order")
pdf("r-test.pdf")
p1
dev.off()
#ggsave("r-test.pdf", plot = p1, width = 20, height = 20, units = "cm")
