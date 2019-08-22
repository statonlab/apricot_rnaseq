setwd("~/Desktop/Jiali/UTK/Peach/")
library(DESeq2)
library("VennDiagram")
library("RColorBrewer")
library("gplots")
library(ggplot2)
library("gridExtra")
library("cowplot")
library(plyr)
library(reshape2)


# Plot mapping metrics
metrics <- read.csv("mapping metrics.csv", header = T)[1:69,]
head(metrics)
metrics$rRNA <- as.numeric(gsub("%","",metrics$rRNA))
metrics <- metrics[order(metrics$rRNA),]
metrics$Sample <- factor(metrics$Sample, levels = metrics$Sample[order(metrics$rRNA)])
metrics$Mapped.reads <- gsub(",","",metrics$Mapped.reads)
metrics$Mapped.reads <- as.numeric(metrics$Mapped.reads)/100000
ggplot(data = metrics) +
  geom_bar(aes(x=Sample,y=Mapped.reads),stat = "identity") +
  geom_point( aes(x=Sample,y=rRNA))+
  scale_y_continuous("Mapped reads", sec.axis=sec_axis(~.*1, name="rRNA %") ) 
#ylab("rRNA (%)")
#  scale_x_discrete(limits=metrics$Sample)

library(latticeExtra)
## 1=== With xyplot, you can easily show both var in the same time :
xyplot(Mapped.reads + rRNA ~ Sample, metrics, type = "l")

## 2=== But it could be nice to have TWO Y axis!

# --> construct separate plots for each series
obj1 <- dotplot(Mapped.reads ~ Sample, metrics , lwd=2)
obj2 <- xyplot(rRNA ~ Sample, metrics, type = "l", lwd=2)

# --> Make the plot with second y axis:
doubleYScale(obj1, obj2, add.ylab2 = TRUE)

## 3=== Same graph with a key legend
png(filename="mapping.png")
doubleYScale(obj1, obj2, text = c("Number of Mapped reads", "Percentage of rRNA") , add.ylab2 = TRUE)
dev.off()


# Read data in, create DEseq dataset.
#-------------------------------------------
mydata <- read.table("gene_count2.0.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
MydataDesign= read.csv("Peach RNA Samples.csv", check.names = F, stringsAsFactors = T, row.names = 1)
colnames(mydata) <- gsub(".Aligned","",colnames(mydata))

# mydata colnames and mydatadesign colnames should be corresponding. Reorder names:
mydata <- mydata[,sort(colnames(mydata))]
colnames(mydata)
MydataDesign <- MydataDesign[sort(row.names(MydataDesign)),]
row.names(MydataDesign)

dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~ Phenotype)
dds <- DESeq(dds)
dds = dds[rowSums(counts(dds))>1, ]
dds <- DESeq(dds)

# PCA plot
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("Timepoint","Genotype"), returnData=TRUE)
pcaData$Timepoint <- gsub("preChill","0", pcaData$Timepoint)
pcaData$Timepoint <- factor(pcaData$Timepoint, c("0","100","600","preBloom-early","1000","preBloom-late"))
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(subset(pcaData, Genotype %in% c("A209","A340","A318","A323")), aes(PC1, PC2, color=Timepoint, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_discrete(breaks=c("A209","A340","A318","A323"))+
  scale_color_discrete(breaks=c("0","100","600","preBloom-early","1000","preBloom-late"))+
  theme_classic(base_size = 20)
#scale_color_manual(values=c("#0033FF", "#6666FF", "#33FF99","#99CC33","#FF9900","#FF3333","#CC33CC"))


# Plot DAM genes
TimePlot <- function(x) {
  data <- plotCounts(dds, gene=x,intgroup=c("Timepoint","Genotype"), returnData=TRUE)
  plot <-
    ggplot(data, aes(x=Timepoint, y=count, color=Genotype, group=Genotype)) +
    geom_point(size=2) + geom_smooth(se = FALSE, method = "loess") +
    scale_color_manual(values=c("#FF3333","#0033FF","#99CCFF","#FF9999")) + scale_y_log10()+
    scale_x_discrete(limits=c("preChill","100","600","preBloom-early","1000","preBloom-late"))
  #print(plot)
  #  return()
}


DAM1 = TimePlot("ppa018667m.g.v1.0")
DAM2 = TimePlot("ppb017585m.g.v1.0")
DAM3 = TimePlot("ppa010758m.g.v1.0")
DAM4 = TimePlot("ppa011123m.g.v1.0")
DAM5 = TimePlot("ppa010822m.g.v1.0")
DAM6 = TimePlot("ppa010714m.g.v1.0")
plot_DAM <- plot_grid(DAM1,DAM2, DAM3, DAM3,DAM5, DAM6, ncol = 3, nrow = 2) + theme(text = element_text(size=10))
#plot_DAM + theme(text = element_text(size=8))
ggsave(filename = "DAM genes.pdf", plot = plot_DAM, width=9, height=3, units="in", scale=3)


# Differential expression test
MydataDesign$phenotype_time <- paste0(MydataDesign$Phenotype,"_",MydataDesign$Timepoint)
head(MydataDesign)
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~ phenotype_time)
dds <- DESeq(dds)
dds = dds[rowSums(counts(dds))>1, ]
dds <- DESeq(dds)
res = results(dds, contrast = c('phenotype_time', 'LowChill_600', 'HighChill_600'), lfcThreshold = 1, alpha = 0.05)
res = results(dds, contrast = c('phenotype_time', 'LowChill_600', 'HighChill_600'))
summary(res)
write.csv(res,'DE_results.csv')

res0 = results(dds, contrast = c("phenotype_time","LowChill_preChill", "HighChill_preChill"), lfcThreshold = 1, alpha = 0.05)
summary(res0)
res100 = results(dds, contrast = c("phenotype_time","LowChill_100", "HighChill_100"), lfcThreshold = 1, alpha = 0.05)
summary(res100)
resBloom = results(dds, contrast = c("phenotype_time","LowChill_preBloom-early", "HighChill_preBloom-late"), lfcThreshold = 1, alpha = 0.05)
summary(resBloom)
res0v600 = results(dds, contrast = c("phenotype_time","LowChill_600", "LowChill_preChill"), lfcThreshold = 1, alpha = 0.05)
summary(res0v600)
res100v600 = results(dds, contrast = c("phenotype_time","LowChill_600", "LowChill_100"), lfcThreshold = 1, alpha = 0.05)
summary(res100v600)

ddsStage = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~ Stage)
ddsStage <- DESeq(ddsStage)
resEndovEco <- results(ddsStage, contrast = c("Stage", "Eco", "Endo"), lfcThreshold = 1, alpha = 0.05)
summary(resEndovEco)

x <- data.frame("Stage" = c("0", "100", "600", "preBloom"), "up" = c(9, 20, 1603, 275), "down" = c(40, 57, 499, 308))
x_melt <- reshape2::melt(x)

ggplot(x_melt, aes(x = Stage, y=value, fill = variable)) +
  geom_bar(stat = "identity") + theme_classic(base_size = 20)+
  scale_fill_manual(values=c("#FF3333","#0033FF")) +
  labs(fill="", x = "Stage", y="Number of genes")


# Venn diagram of apricot and peach data
res600 <- res[ which(res$padj <= 0.05),]
res_1v2 <- read.csv("../Apricot/Sig_genes_time_1v2_new.csv",row.names = 1,header = T)

v_1 <- row.names(res600)
v_2 <- row.names(res_1v2)

# venn diagram - ALL
area1=length(v_1)
area2=length(v_2)
#---pairs
n12=length(intersect(v_1,v_2))
#-----------------------------------
grid.newpage()
draw.pairwise.venn(area1, area2, n12,
                   scaled=T, 
                   fill=c("red", "yellow"))

# Generate files for GO enrichment
library(plyr)
library(tidyr)

annotat <- read.csv("../Apricot/Ppersica_298_v2.1.annotation_info.txt", header = TRUE, comment.char = '', sep='\t')
write.csv(x = res600, file = "Sig_genes_600.csv")

res_Sig600 <- read.csv("Sig_genes_600.csv",header = TRUE)
res_Sig600$X = gsub('\\.v2\\.1', '', res_Sig600$X)
colnames(res_Sig600)[1] <- "locusName"
match1 <- join(res_Sig600, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
locusname <- read.csv("../Apricot/Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_600.csv")


Gene_up <- match3[which(match3$log2FoldChange>0),]
Gene_down <- match3[which(match3$log2FoldChange<0),]
#head(match3$GO)
Gene_up$GO <- as.character(Gene_up$GO)
s <- strsplit(Gene_up$GO, split = ",")
GO_600up <- data.frame(locus = rep(Gene_up$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_600up, file = "GO_1v2up", sep = " ", quote = F, row.names = F, col.names = F)

Gene_down$GO <- as.character(Gene_down$GO)
s <- strsplit(Gene_down$GO, split = ",")
GO_600down <- data.frame(locus = rep(Gene_down$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_600down, file = "GO_1v2down", sep = " ", quote = F, row.names = F, col.names = F)


# overlap genes GO term
inter <- res600[which(rownames(res600) %in% intersect(v_1,v_2)),]
write.csv(x = inter, file = "Sig_genes_inter.csv")

res_Siginter <- read.csv("Sig_genes_inter.csv",header = TRUE)
res_Siginter$X = gsub('\\.v2\\.1', '', res_Siginter$X)
colnames(res_Siginter)[1] <- "locusName"
match1 <- join(res_Siginter, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
locusname <- read.csv("../Apricot/Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
write.csv(match3,file = "DEgene_inter.csv")

match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_inter <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_inter, file = "GO_inter", sep = " ", quote = F, row.names = F, col.names = F)


## Add rRNA as a factor
# create data discription
dist_rld_t.env <- MydataDesign
dist_rld_t.env$Sample <- rownames(MydataDesign)
dist_rld_t.env <- merge(dist_rld_t.env, metrics, by = "Sample")
dist_rld_t.env$rRNA <- dist_rld_t.env$rRNA/100
rownames(dist_rld_t.env) <- dist_rld_t.env$Sample
dist_rld_t.env <- dist_rld_t.env[,-1]

MydataDesign_rRNA <- dist_rld_t.env[,-c(5)]
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_rRNA, design = ~ Phenotype)
dds = dds[rowSums(counts(dds))>1, ]
dds <- DESeq(dds)

# PCA plot, rRNA content as a transparency factor
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = c("Timepoint","Genotype","rRNA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(subset(pcaData, Genotype %in% c("A209","A340","A318","A323")), aes(PC1, PC2, color=Timepoint, shape=Genotype)) +
  geom_point(size=3, aes(alpha=1-rRNA)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_shape_discrete(breaks=c("A209","A340","A318","A323"))+
  scale_color_discrete(breaks=c("preChill","100","600","preBloom-early","1000","preBloom-late"))+
  theme(text = element_text(size=16))
#scale_color_manual(values=c("#0033FF", "#6666FF", "#33FF99","#99CC33","#FF9900","#FF3333","#CC33CC"))

library('Rmisc')
data <- plotCounts(dds, gene="Prupe.3G128200.v2.1",intgroup=c("Timepoint","Genotype"), returnData=TRUE)
datasum <- summarySE(data, measurevar = "count", groupvars = c("Timepoint","Genotype"))
plot <-
  ggplot(datasum, aes(x=Timepoint, y=count, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=count-se,ymax=count+se), width=0.1)+
  scale_color_manual(values=c("#FF3333","#0033FF","#99CCFF","#FF9999")) + #scale_y_log10()+
  scale_x_discrete(limits=c("preChill","100","600","preBloom-early","1000","preBloom-late"))
plot
