#!/Users/sean/anaconda3/envs/renv/bin/Rscript

# #/usr/bin/env Rscript


library(ggpubr)
library(ggplot2)

dat1<-read.csv("R-input-data.tsv", sep="\t", header=T)

set_shape_to_genometype <- function(dat1) {
    if ("Query" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Query") == 1) {
        sympos1 <- 24
    }}
    if ("Ref-Isolate" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-Isolate") == 1) {
        sympos1 <- 22
    }}
    if ("Ref-MAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-MAG") == 1) {
        sympos1 <- 21
    }}
    if ("Ref-SAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-SAG") == 1) {
        sympos1 <- 23
    }}
    if ("Ref-Isolate" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-Isolate") == 2) {
        sympos2 <- 22
    }}
    if ("Ref-MAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-MAG") == 2) {
        sympos2 <- 21
    }}
    if ("Ref-SAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-SAG") == 2) {
        sympos2 <- 23
    }}
    if ("Ref-MAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-MAG") == 3) {
        sympos3 <- 21
    }}
    if ("Ref-SAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-SAG") == 3) {
        sympos3 <- 23
    }}
    if ("Ref-SAG" %in% levels(dat1$GenomeType)) {
    if (which(levels(dat1$GenomeType) == "Ref-SAG") == 4) {
        sympos4 <- 23
    }}
    if (length(levels(dat1$GenomeType)) == 4) {
        plotsym <- scale_shape_manual(values=c(sympos1, sympos2, sympos3, sympos4))
    }
    if (length(levels(dat1$GenomeType)) == 3) {
        plotsym <- scale_shape_manual(values=c(sympos1, sympos2, sympos3))
    }
    if (length(levels(dat1$GenomeType)) == 2) {
        plotsym <- scale_shape_manual(values=c(sympos1, sympos2))
    }
    if (length(levels(dat1$GenomeType)) == 1) {
    	   print("1 detected")
        plotsym <- scale_shape_manual(values=c(sympos1))
    }
    return(plotsym)
}

set_shape_to_genometype(dat1)

p1 <- ggplot(dat1, aes(x=PC1, y=PC2))
p1 <- p1 + theme_bw()
p1 <- p1 + geom_point(aes(fill=Class, shape=GenomeType))
p1 <- p1 + set_shape_to_genometype(dat1)
p1 <- p1 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p1 <- p1 + guides(shape = guide_legend(title="Genome Type"))
p1 <- p1 + xlab("PC1") + ylab("PC2")
p1 <- p1 + labs(title = "PCA of RAST annotation", subtitle = "Color: Class taxonomy level; Shape: Genome Type")

p2 <- ggplot(dat1, aes(x=PC1, y=PC2))
p2 <- p2 + theme_bw()
p2 <- p2 + geom_point(aes(fill=Order, shape=GenomeType))
p2 <- p2 + set_shape_to_genometype(dat1)
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
p5 <- p5 + set_shape_to_genometype(dat1)
p5 <- p5 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p5 <- p5 + guides(shape = guide_legend(title="Genome Type"))
p5 <- p5 + xlab("Completeness (%)") + ylab("Contamination (%)")
p5 <- p5 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Class taxonomy level; Shape: Genome Type")


p6 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p6 <- p6 + theme_bw()
p6 <- p6 + geom_point(aes(fill=Order, shape=GenomeType))
p6 <- p6 + set_shape_to_genometype(dat1)
p6 <- p6 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p6 <- p6 + guides(shape = guide_legend(title="Genome Type"))
p6 <- p6 + xlab("Completeness (%)") + ylab("Contamination (%)")
p6 <- p6 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")


p7 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p7 <- p7 + theme_bw()
p7 <- p7 + geom_point(aes(fill=Class, shape=GenomeType))
p7 <- p7 + set_shape_to_genometype(dat1)
p7 <- p7 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p7 <- p7 + guides(shape = guide_legend(title="Genome Type"))
p7 <- p7 + xlab("Completeness (%)") + ylab("Contamination (%)")
p7 <- p7 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Class taxonomy level; Shape: Genome Type")
p7 <- p7 + facet_wrap(~Class, nrow = length(levels(dat1$Class)))
p7 <- p7 + ylim(0,20) + xlim(0,100)
if (length(levels(dat1$Class)) > 10) {
    p7 <- p7 + theme(strip.text.x = element_blank())
}



p8 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p8 <- p8 + theme_bw()
p8 <- p8 + geom_point(aes(fill=Order, shape=GenomeType))
p8 <- p8 + set_shape_to_genometype(dat1)
p8 <- p8 + guides(fill = guide_legend(title="Class", override.aes=list(shape=21)))
p8 <- p8 + guides(shape = guide_legend(title="Genome Type"))
p8 <- p8 + xlab("Completeness (%)") + ylab("Contamination (%)")
p8 <- p8 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")
p8 <- p8 + facet_wrap(~Order, nrow = length(levels(dat1$Order)))
p8 <- p8 + ylim(0,20) + xlim(0,100)
if (length(levels(dat1$Order)) > 10) {
    p8 <- p8 + theme(strip.text.x = element_blank())
}


p9 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p9 <- p9 + theme_bw()
p9 <- p9 + geom_rect(aes(xmin=0, xmax = 50, ymin = 0, ymax = 2), fill = "pink", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=0, xmax = 50, ymin = 2, ymax = Inf), fill = "red", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=50, xmax = 75, ymin = 4, ymax = Inf), fill = "red", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=50, xmax = 75, ymin = 2, ymax = 4), fill = "pink", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=50, xmax = 75, ymin = 0, ymax = 2), fill = "yellow", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=75, xmax = 100, ymin = 4, ymax = Inf), fill = "pink", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=75, xmax = 87.5, ymin = 2, ymax = Inf), fill = "pink", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=87.5, xmax = 100, ymin = 2, ymax = 4), fill = "yellow", alpha = 0.01)
p9 <- p9 + geom_rect(aes(xmin=75, xmax = 100, ymin = 0, ymax = 2), fill = "green", alpha = 0.01)
p9 <- p9 + geom_point(aes(fill=Order, shape=GenomeType))
p9 <- p9 + set_shape_to_genometype(dat1)
p9 <- p9 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p9 <- p9 + guides(shape = guide_legend(title="Genome Type"))
p9 <- p9 + xlab("Completeness (%)") + ylab("Contamination (%)")
p9 <- p9 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")
p9 <- p9 + ylim(0,20) + xlim(0,100)


p10 <- ggplot(dat1, aes(x=Completeness, y=Contamination))
p10 <- p10 + theme_bw()
p10 <- p10 + geom_rect(aes(xmin=0, xmax = 50, ymin = 0, ymax = 2), fill = "pink", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=0, xmax = 50, ymin = 2, ymax = Inf), fill = "red", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=50, xmax = 75, ymin = 4, ymax = Inf), fill = "red", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=50, xmax = 75, ymin = 2, ymax = 4), fill = "pink", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=50, xmax = 75, ymin = 0, ymax = 2), fill = "yellow", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=75, xmax = 100, ymin = 4, ymax = Inf), fill = "pink", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=75, xmax = 87.5, ymin = 2, ymax = Inf), fill = "pink", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=87.5, xmax = 100, ymin = 2, ymax = 4), fill = "yellow", alpha = 0.01)
p10 <- p10 + geom_rect(aes(xmin=75, xmax = 100, ymin = 0, ymax = 2), fill = "green", alpha = 0.01)
p10 <- p10 + geom_point(aes(fill=Order, shape=GenomeType))
p10 <- p10 + set_shape_to_genometype(dat1)
p10 <- p10 + guides(fill = guide_legend(title="Order", override.aes=list(shape=21)))
p10 <- p10 + guides(shape = guide_legend(title="Genome Type"))
p10 <- p10 + xlab("Completeness (%)") + ylab("Contamination (%)")
p10 <- p10 + labs(title = "Genome Completeness and Contamination Scatter Plot", subtitle = "Color: Order taxonomy level; Shape: Genome Type")
p10 <- p10 + ylim(0,20) + xlim(0,100)



p_final1 <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2)
p_final2 <- ggarrange(p3, p4, labels = c("C", "D"), ncol = 1, nrow = 2)
p_final3 <- ggarrange(p5, p6, labels = c("E", "F"), ncol = 1, nrow = 2)
p_final4 <- ggarrange(p7, labels = c("G"), ncol = 1, nrow = 1)
p_final5 <- ggarrange(p8, labels = c("H"), ncol = 1, nrow = 1)
p_final6 <- ggarrange(p9, p10, labels = c("I", "J"), ncol = 1, nrow = 2)

pdf("R-output-plot.pdf", width=8.5,height=11)
p_final1
p_final2
p_final3
p_final4
p_final5
p_final6
dev.off()
