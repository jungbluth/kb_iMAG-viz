#!/Users/sean/anaconda3/envs/renv/bin/Rscript

library(ggplot2)
library(ggpubr)

dat1<-read.csv("R-input-data.tsv", sep="\t", header=T)


p1 <- ggplot(dat1, aes(x=PC1, y=PC2))
p1 <- p1 + theme_bw()
p1 <- p1 + geom_point(aes(fill=Class, shape=GenomeType))
p1 <- p1 + scale_shape_manual(values=c(24, 22, 21, 23))
p1 <- p1 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p1 <- p1 + guides(shape = guide_legend(title="Genome Type"))
p1 <- p1 + xlab("PC1") + ylab("PC2")
p1 <- p1 + labs(title = "PCA of RAST annotation", subtitle = "Color: Class taxonomy level; Shape: Genome Type")

p2 <- ggplot(dat1, aes(x=PC1, y=PC2))
p2 <- p2 + theme_bw()
p2 <- p2 + geom_point(aes(fill=Order, shape=GenomeType))
p2 <- p2 + scale_shape_manual(values=c(24, 22, 21, 23))
p2 <- p2 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p2 <- p2 + guides(shape = guide_legend(title="Genome Type"))
p2 <- p2 + xlab("PC1") + ylab("PC2")
p2 <- p2 + labs(title = "PCA of RAST annotation", subtitle = "Color: Order taxonomy level; Shape: Genome Type")

p3 <- ggplot(dat1, aes(x=PC1, y=PC2))
p3 <- p3 + theme_bw()
p3 <- p3 + geom_point(aes(color=Class, size=Completeness))
p3 <- p3 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p3 <- p3 + guides(size = guide_legend(title="Completeness (%)"))
p3 <- p3 + xlab("PC1") + ylab("PC2")
p3 <- p3 + labs(title = "PCA of RAST annotation", subtitle = "Color by Class taxonomy level; Size by Completeness")


p4 <- ggplot(dat1, aes(x=PC1, y=PC2))
p4 <- p4 + theme_bw()
p4 <- p4 + geom_point(aes(color=Class, size=Contamination))
p4 <- p4 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p4 <- p4 + guides(size = guide_legend(title="Contamination (%)"))
p4 <- p4 + xlab("PC1") + ylab("PC2")
p4 <- p4 + labs(title = "PCA of RAST annotation", subtitle = "Color by Order taxonomy level; Size by Contamination")


p5 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p5 <- p5 + theme_bw()
p5 <- p5 + geom_point(aes(fill=Class, shape=GenomeType))
p5 <- p5 + scale_shape_manual(values=c(24, 22, 21, 23))
p5 <- p5 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p5 <- p5 + guides(shape = guide_legend(title="Genome Type"))
p5 <- p5 + xlab("Completeness (%)") + ylab("Contamination (%)")
p5 <- p5 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Class taxonomy level; Shape: Genome Type")


p6 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p6 <- p6 + theme_bw()
p6 <- p6 + geom_point(aes(fill=Order, shape=GenomeType))
p6 <- p6 + scale_shape_manual(values=c(24, 22, 21, 23))
p6 <- p6 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p6 <- p6 + guides(shape = guide_legend(title="Genome Type"))
p6 <- p6 + xlab("Completeness (%)") + ylab("Contamination (%)")
p6 <- p6 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")


p7 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p7 <- p7 + theme_bw()
p7 <- p7 + geom_rect(aes(xmin=0, xmax = 50, ymin = 0, ymax = 2), fill = "pink", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=0, xmax = 50, ymin = 2, ymax = Inf), fill = "red", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=50, xmax = 75, ymin = 4, ymax = Inf), fill = "red", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=50, xmax = 75, ymin = 2, ymax = 4), fill = "pink", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=50, xmax = 75, ymin = 0, ymax = 2), fill = "yellow", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=75, xmax = 100, ymin = 4, ymax = Inf), fill = "pink", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=75, xmax = 87.5, ymin = 2, ymax = Inf), fill = "pink", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=87.5, xmax = 100, ymin = 2, ymax = 4), fill = "yellow", alpha = 0.01)
p7 <- p7 + geom_rect(aes(xmin=75, xmax = 100, ymin = 0, ymax = 2), fill = "green", alpha = 0.01)
p7 <- p7 + geom_point(aes(fill=Order, shape=GenomeType))
p7 <- p7 + scale_shape_manual(values=c(24, 22, 21, 23))
p7 <- p7 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p7 <- p7 + guides(shape = guide_legend(title="Genome Type"))
p7 <- p7 + xlab("Completeness (%)") + ylab("Contamination (%)")
p7 <- p7 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")


p8 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p8 <- p8 + theme_bw()
p8 <- p8 + geom_rect(aes(xmin=0, xmax = 50, ymin = 0, ymax = 2), fill = "pink", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=0, xmax = 50, ymin = 2, ymax = Inf), fill = "red", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=50, xmax = 75, ymin = 4, ymax = Inf), fill = "red", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=50, xmax = 75, ymin = 2, ymax = 4), fill = "pink", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=50, xmax = 75, ymin = 0, ymax = 2), fill = "yellow", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=75, xmax = 100, ymin = 4, ymax = Inf), fill = "pink", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=75, xmax = 87.5, ymin = 2, ymax = Inf), fill = "pink", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=87.5, xmax = 100, ymin = 2, ymax = 4), fill = "yellow", alpha = 0.01)
p8 <- p8 + geom_rect(aes(xmin=75, xmax = 100, ymin = 0, ymax = 2), fill = "green", alpha = 0.01)
p8 <- p8 + geom_point(aes(fill=Order, shape=GenomeType))
p8 <- p8 + scale_shape_manual(values=c(24, 22, 21, 23))
p8 <- p8 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p8 <- p8 + guides(shape = guide_legend(title="Genome Type"))
p8 <- p8 + xlab("Completeness (%)") + ylab("Contamination (%)")
p8 <- p8 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")


p_final1 <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2)
p_final2 <- ggarrange(p3, p4, labels = c("C", "D"), ncol = 1, nrow = 2)
p_final3 <- ggarrange(p5, p6, labels = c("E", "F"), ncol = 1, nrow = 2)
p_final4 <- ggarrange(p7, p8, labels = c("G", "H"), ncol = 1, nrow = 2)

pdf("R-output-plot.pdf", width=8.5,height=11)
p_final1
p_final2
p_final3
p_final4
dev.off()
