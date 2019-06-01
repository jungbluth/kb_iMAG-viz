#!/usr/bin/env Rscript

library(ggplot2)
dat1<-read.csv("_Plot-input-data.tsv", sep="\t", header=T)
head(dat1)
p1 <- ggplot(dat1) + theme_bw() + geom_point(aes(x=dat1$PC1,y=dat1$PC2,alpha=0.8))

pdf("r-test.pdf")
p1
dev.off()
#ggsave("r-test.pdf", plot = p1, width = 20, height = 20, units = "cm")
