#!/usr/bin/env Rscript

library(ggplot2)
dat1<-read.csv("R-input-data.tsv", sep="\t", header=T)
head(dat1)
p1 <- ggplot(dat1, aes(x=PC1, y=PC2))
p1 <- p1 + theme_bw()
p1 <- p1 + geom_point(aes(fill=Class, shape=GenomeType))
p1 <- p1 + scale_shape_manual(values=c(24, 22, 21, 23))
p1 <- p1 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p1 <- p1 + xlab("PC1") + ylab("PC2")
pdf("R-output-plot.pdf")
p1
dev.off()
