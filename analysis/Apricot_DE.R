library(DESeq2)
install.packages('VennDiagram')
library("VennDiagram")
library("RColorBrewer")
library("gplots")
library(ggplot2)
install.packages("gridExtra")
library("gridExtra")
install.packages("cowplot")
library("cowplot")
library(plyr)
library(reshape2)

# Load HTseq count data
mydata <- read.table("gene_counts.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
head(mydata)

# Import data description table, regroup 800 timepoint, early bloom tree as flower-800, late bloom tree as bud-800
MydataDesign_1 <- read.csv("MydataDesign_1.csv",check.names = F, stringsAsFactors = T, row.names = 1)
MydataDesign_1$phenotype_stage <- factor(paste0(MydataDesign_1$phenotype,"_",MydataDesign_1$Stage))
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_1, design = ~Stage)
dim(dds)
dds = DESeq(dds)

# PCA plot
rld <- rlog(dds, blind = FALSE)
#pcahead(assay(rld), 3)
pcaData <- plotPCA(rld, intgroup = c("Stage", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(subset(pcaData, genotype %in% c("A1956","A2137","A1267","A660")), aes(PC1, PC2, color=Stage, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(text = element_text(size=16))+
  theme_classic()+
  scale_shape_manual(values = c(16,4,15,8)) +
  scale_color_manual(values=c("#0033FF", "#6666FF", "#33FF99","#99CC33","#FF9900","#FF3333","#CC33CC"))
ggsave("PCAplot0423.png")  

ggplot(pcaData, aes(PC1, PC2, color=timepoint)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  facet_wrap(~phenotype)
ggplot(pcaData, aes(PC1, PC2)) +
  geom_line() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  facet_wrap(~phenotype) +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggplot(pcaDataGeno, aes(PC1, PC2)) +
  geom_line() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  facet_wrap(~genotype) +
  theme(strip.background = element_blank(), strip.placement = "outside")

# a new design based on PCA clusters group into 1_endo, 2_eco, 3_sepal, 4_petal, 5_flower
MydataDesign_new = read.csv("Mydatadesign.csv", header = T, row.names = 1)
MydataDesign_new$timepoint <- factor(MydataDesign_new$timepoint)
MydataDesign_new$time_genotype <- factor(paste0(MydataDesign_new$timepoint,'_',MydataDesign_new$genotype))
MydataDesign_new$time_phenotype <- factor(paste0(MydataDesign_new$timepoint,'_',MydataDesign_new$phenotype))
head(MydataDesign_new)

ddsEndo = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_new, design = ~ time_phenotype)
ddsEndo = DESeq(ddsEndo)
# test the genotypic difference among the endo
resTime_endo <- results(ddsEndo, contrast = c("time_phenotype", "1_early","1_late"), alpha = 0.05, lfcThreshold = 1)
summary(resTime_endo)

# Test the DE comparing endo and eco samples based on PCA
ddsPhase = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_new, design = ~ timepoint)
ddsPhase <- DESeq(ddsPhase)
resTime_1v2 <- results(ddsPhase, contrast = c("timepoint", "1","2"), alpha = 0.05, lfcThreshold = 1)
summary(resTime_1v2)

write.csv(x = resTime_1v2, file = "Sig_genes_time_1v2_new.csv")
resTime_1v2 <- read.csv("Sig_genes_time_1v2_new.csv",header = TRUE)
colnames(resTime_1v2)[1] <- "locusName"
resTime_1v2$locusName = gsub('\\.v2\\.1', '', resTime_1v2$locusName)
match1 <- join(resTime_1v2, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
write.csv(match3, file = "DEgene_time1v2.csv")


## Plot flower locus T gene
data <- plotCounts(dds, gene="Prupe.1G398700.v2.1",intgroup=c("Stage","genotype"), returnData=TRUE)
plot <- ggplot(data, aes(x=Stage, y=count, color=genotype, group=genotype)) +
  geom_point(size=2) + geom_smooth(se = FALSE, method = "loess") #+ scale_y_log10() # geom_line(size=1)
plot + theme(text = element_text(size=10), 
               plot.margin = unit(c(2, 2, 2, 2), "cm")) + 
  labs(title = "FLOWER LOCUS T") + scale_color_manual(values=c("#0033FF","#FF3333","#FF9999","#99CCFF"))

#-------------------------
#######################################################################################
####Compare DE gene at each time point and stage
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_1, design = ~phenotype_stage)
dds = ddsC[rowSums(counts(dds))>1, ]
dim(dds)
dds = DESeq(dds)

# 0 chill hour
res0 <- results(dds, contrast = c('phenotype_stage', 'early_0', 'late_0'),lfcThreshold = 1, alpha = 0.05)
summary(res0)
plotMA(res0, main = "Early vs late at 0 chill hour")

res0Sig <- res0[ which(res0$padj < 0.05),]
res0Sig
write.csv(x = res0Sig, file = "Sig_genes_0.csv")
res0Sig <- read.csv("Sig_genes_0.csv",header = TRUE)

res0Sig$X = gsub('\\.v2\\.1', '', res0Sig$X)
colnames(res0Sig)[1] <- "locusName"
match1 <- join(res0Sig, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_0.csv")

head(match3$GO)
match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_match3 <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_match3, file = "GO_DEgene0", sep = " ", quote = F, row.names = F, col.names = F)

####################################################
# 100 chill hour
res100 <- results(dds, contrast = c('phenotype_stage', 'early_100', 'late_100'), lfcThreshold = 1, alpha = 0.05)
summary(res100)
plotMA(res100, main = "Early flowering vs Late flowering at 100 chill hours")

res100Sig <- res100[ which(res100$padj < 0.05),]
write.csv(x = res100Sig, file = "Sig_genes_100.csv")
res100Sig <- read.csv("Sig_genes_100.csv",header = TRUE)

res100Sig$X = gsub('\\.v2\\.1', '', res100Sig$X)
colnames(res100Sig)[1] <- "locusName"
match1 <- join(res100Sig, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_100.csv")

head(match3$GO)
match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_match3 <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_match3, file = "GO_DEgene100", sep = " ", quote = F, row.names = F, col.names = F)

#########################################################
# 400 chill hour
res400 <- results(dds, contrast = c('phenotype_stage', 'early_400', 'late_400'), lfcThreshold = 1, alpha = 0.05)
summary(res400)
plotMA(res400, main = "Early flowering vs Late flowering at 400 chill hour") 

res400Sig <- res400[ which(res400$padj < 0.05),]
res400Sig
write.csv(x = res400Sig, file = "Sig_genes_400.csv")
res400Sig <- read.csv("Sig_genes_400.csv",header = TRUE)
res400Sig$X = gsub('\\.v2\\.1', '', res400Sig$X)
colnames(res400Sig)[1] <- "locusName"
match1 <- join(res400Sig, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_400.csv")

head(match3$GO)
match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_match3 <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_match3, file = "GO_DEgene400", sep = " ", quote = F, row.names = F, col.names = F)

# compare 0 vs 100 chill hours
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign_1, design = ~Stage)
dim(dds)
dds = DESeq(dds)

res0v100 <- results(dds, contrast = c('Stage', '0', '100'), lfcThreshold = 1, alpha = 0.05)
summary(res0v100)
plotMA(res0v100, main = "0 vs 100 chill hour") 

res0v100Sig <- res0v100[ which(res0v100$padj < 0.05),]
res0v100Sig
write.csv(x = res0v100Sig, file = "Sig_genes_0v100.csv")
res0v100Sig <- read.csv("Sig_genes_0v100.csv",header = TRUE)
res0v100Sig$X = gsub('\\.v2\\.1', '', res0v100Sig$X)
colnames(res0v100Sig)[1] <- "locusName"
match1 <- join(res0v100Sig, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_0v100.csv")

head(match3$GO)
match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_match3 <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_match3, file = "GO_DEgene0v100", sep = " ", quote = F, row.names = F, col.names = F)

# compare 100 vs 400 chill hours
res100v400 <- results(dds, contrast = c('Stage', '100', '400'), lfcThreshold = 1, alpha = 0.05)
summary(res100v400)
plotMA(res100v400, main = "100 vs 400 chill hour") 

res100v400Sig <- res400[ which(res100v400$padj < 0.05),]
res100v400Sig
write.csv(x = res100v400Sig, file = "Sig_genes_100v400.csv")
res100v400Sig <- read.csv("Sig_genes_100v400.csv",header = TRUE)
res100v400Sig$X = gsub('\\.v2\\.1', '', res100v400Sig$X)
colnames(res100v400Sig)[1] <- "locusName"
match1 <- join(res100v400Sig, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
write.csv(match3, file = "DEgene_100v400.csv")

head(match3$GO)
match3$GO <- as.character(match3$GO)
s <- strsplit(match3$GO, split = ",")
GO_match3 <- data.frame(locus = rep(match3$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_match3, file = "GO_DEgene100v400", sep = " ", quote = F, row.names = F, col.names = F)

#========================================================
# Gene clustering
library("genefilter")
library( "RColorBrewer" )
library("gplots")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 2000)
#----------------
matricx <- assay(rld)[ topVarGenes, ]
matricx_1 <- matricx[,grep('C|D', colnames(matricx))]
head(matricx_1)

heatmap.2( matricx_1, scale="row",
           trace="none", dendrogram="column", labRow = F,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
#--------------
TopEcoVarGenes <- head(order(Time_1v2$padj, decreasing = T),500)
Variables <- rownames(MydataDesign_new[MydataDesign_new$timepoint==1 | MydataDesign_new$timepoint==2,])
heatmap.2( assay(rld)[ TopEcoVarGenes, c(Variables)], scale="row",
           trace="none", dendrogram="column", labRow = F,
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
label <- as.data.frame(colData(rld)[, c("timepoint", "phenotype")])


## ------------Analysis 05232018----------------------------
# Venn diagram
DE0 <- read.csv("DEgene_0.csv", header = T)
DE100 <- read.csv("DEgene_100.csv", header = T)
DE400 <- read.csv("DEgene_400.csv", header = T)
DE0v100 <- read.csv("DEgene_0v100.csv", header = T)
DE100v400 <- read.csv("DEgene_100v400.csv", header = T)
DE1v2 <- read.csv("DEgene_time1v2.csv", header = T)

Gene_0 <- DE0$locusName
Gene_100 <- DE100$locusName
Gene_400 <- DE400$locusName
Gene_0v100 <- DE0v100$locusName
Gene_100v400 <- DE100v400$locusName
Gene_1v2 <- DE1v2$locusName

# venn diagram - quatria
area1=length(Gene_0)
area2=length(Gene_100)
area3=length(Gene_400)
area4=length(Gene_1v2)
#---pairs
n12=length(intersect(Gene_0,Gene_100))
n13=length(intersect(Gene_0,Gene_400))
n14=length(intersect(Gene_0,Gene_1v2))
n23=length(intersect(Gene_100,Gene_400))
n24=length(intersect(Gene_100,Gene_1v2))
n34=length(intersect(Gene_400,Gene_1v2))
#---trios
n123=length(Reduce(intersect,list(Gene_0,Gene_100,Gene_400)))
n124=length(Reduce(intersect,list(Gene_0,Gene_0v100,Gene_1v2)))
n134=length(Reduce(intersect,list(Gene_0,Gene_400,Gene_1v2)))
n234=length(Reduce(intersect,list(Gene_100,Gene_400,Gene_1v2)))
#---quats
n1234=length(Reduce(intersect,list(Gene_0,Gene_100,Gene_400,Gene_1v2)))
#-----------------------------------
grid.newpage()
draw.quad.venn(area1, area2, area3, area4, n12, n13, n14, n23, n24,
               n34, n123, n124, n134, n234, n1234)

v1 = draw.pairwise.venn(area1, area4, n14,
                        scaled=FALSE, 
                        fill=c("red", "yellow"),
                        c("1v2","0"))
grid.newpage()
v2 = draw.pairwise.venn(area2, area4, n24,
                        scaled=FALSE, 
                        fill=c("red", "yellow"),
                        c("1v2","100"))
grid.newpage()
v3 = draw.pairwise.venn(area3, area4, n34,
                        scaled=FALSE, 
                        fill=c("red", "yellow"),
                        c("1v2","400"))

## ---0v100, 100v400
area1=length(Gene_0v100)
area2=length(Gene_100v400)
area3=length(Gene_1v2)

#---pairs
n12=length(intersect(Gene_0v100,Gene_100v400))
n13=length(intersect(Gene_0v100,Gene_1v2))
n23=length(intersect(Gene_100v400,Gene_1v2))

#---trios
n123=length(Reduce(intersect,list(Gene_0v100,Gene_100v400,Gene_1v2)))

grid.newpage()
draw.triple.venn(area1, area2, area3, n12, n13, n23, n123,
                 scaled=FALSE, 
                 fill=c("red", "blue","yellow"), 
                 c("0v100", "100v400","1v2"))
