setwd("~/Desktop/Jiali/UTK/Peach/")
library(DESeq2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("WGCNA") # install, ERROR: compilation failed for package ‘WGCNA’, download the source file to install
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")) 
library(WGCNA)
options(stringsAsFactors = F)
# Normalize peach and apricot datasets
PeachData = read.table("gene_count2.0.txt", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
ApricotData <- read.table("../Apricot/gene_counts.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
colnames(PeachData) <- gsub(".Aligned","",colnames(PeachData))
PeachData <- PeachData[,sort(colnames(PeachData))]
colnames(PeachData)
combData <- cbind(PeachData,ApricotData)
MydataDesign <- read.csv("Both Datadesign.csv", header = T, row.names = 1)
dds = DESeqDataSetFromMatrix(countData = combData, colData = MydataDesign, design = ~ Stage)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))
write.table(normalized_counts, "comb normalized counts.txt")

normalized_counts1 <- read.table("comb normalized counts.txt", header = T, row.names = 1)
names(normalized_counts1) <- gsub("X","",names(normalized_counts1))
# remove genes with more than half 0 values
normalized_counts1[normalized_counts1 == 0] <- NA
data <- normalized_counts1[ -which(rowMeans(is.na(normalized_counts1)) > 0.5), ]
data[is.na(data)] <- 0

datExpr0 = as.data.frame(t(data))
datExpr0 <- as.data.frame(lapply(datExpr0,as.numeric)) #datExpr0 has to be numeric.
rownames(datExpr0) <- colnames(combData)
datExpr0[1:3,1:3]
dim(datExpr0)

# check missing value and outliers
gsg = goodSamplesGenes(datExpr0,  verbose = 3)
gsg$allOK
# Returns True means all genes pass the cuts, if not, the following scripts can remove the outliers.
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Auto one step construction
allowWGCNAThreads() # Work in Rstudio

# Choose a set of soft-thresholding powers
powers = c(c(1:16))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# One-step network construction and module detection
net = blockwiseModules(datExpr0, power = 9,
                       TOMType = "unsigned", minModuleSize = 200,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "bothTOM",
                       verbose = 3)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs1 = net$MEs;

datTraits0 = read.csv("bothTrait.csv", header = T,row.names = 1,stringsAsFactors = T)
rownames(datTraits0) 
rownames(datExpr0) <- colnames(data)

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits0[,-12], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(11,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf("corrplot_0428.pdf",width = 6,height = 6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0)[-12],
               yLabels = names(MEs1),
               ySymbols = names(MEs1),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-factors relationships"))
dev.off()

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 9)

## network visualization
annot = read.csv(file = "../Apricot/Ppersica_298_v2.1.annotation_info.txt",header = T,sep = "\t");
annot$arabi.hit = gsub('.{2}$',"",annot$Best.hit.arabi.name)
annot$arabi.hit = gsub("T","t",annot$arabi.hit)
annot$arabi.hit = gsub("G","g",annot$arabi.hit)

## export to cytoscape
# Select modules
modules = "6";
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleLabels, modules));
modProbes = probes[inModule];
modProbes = gsub("\\.v2.1",'',modProbes)
modGenes = annot$arabi.hit[match(modProbes, annot$locusName)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("co_CytoscapeInput-edges-ME5.txt", sep=""),
                               nodeFile = paste("co_CytoscapeInput-nodes-ME5.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# DEgenes from DEseq
DEprof = read.csv("DEgene_inter.csv")
DEgenes = DEprof$locusName
# Genes in each module
ModuleTab <- as.data.frame(matrix(NA,nrow=16,ncol=2))
names(ModuleTab) <- c("Non-DEgenes","DEgenes")
row.names(ModuleTab) <- names(MEs1)
for (i in 0:15) {
  inModule = (moduleLabels==i)
  modProbes = probes[inModule]
  modProbes = gsub("\\.v2.1",'',modProbes)
  Total = length(modProbes)
  interDEgenes = length(intersect(DEgenes,modProbes))
  noDEgenes = Total - interDEgenes
  r= which(substring(names(MEs1),3)==i)
  ModuleTab[r,1] = noDEgenes
  ModuleTab[r,2] = interDEgenes
} 
# Plot bar chart 
tModuleTab <- as.data.frame(t(ModuleTab))
revdat <- tModuleTab[,rev(names(tModuleTab))]
dat = as.matrix(revdat)

sizeGrWindow(11,12)

pdf( "conetwork/MEbar_0430.pdf", width = 12, height = 6 )
par(mfrow = c(1,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0)[-12],
               yLabels = names(MEs1),
               ySymbols = names(MEs1),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-factors relationships"))

barplot(dat, main="Number of genes in each module",
        xlab="Number of genes", las=2, col=c("blue","red"),
        horiz = T, legend = names(ModuleTab),args.legend = list(x="topright"))
dev.off()

# Expression profile of each module
# 1. the genes inside a given module are summarize with the module eigengene, which can be considered as the best summary of the standardized module expression data
MydataDesign <- read.csv("Both Datadesign.csv", header = T, row.names = 1)
datAll <- cbind(MEs1, MydataDesign[,c(1,3)])
datAll$Timepoint <- gsub("C", "Sepal", datAll$Timepoint)
datAll$Timepoint <- gsub("D", "Petal", datAll$Timepoint)
PeachAll <- datAll[1:69,]
ApricotAll <- datAll[70:129,]
PeachAll$Genotype <- factor(PeachAll$Genotype, levels = c("A209", "A340", "A318", "A323"))
PeachAll$Timepoint <- gsub("preChill","0", PeachAll$Timepoint)
datasum <- summarySE(PeachAll, measurevar = "ME11", groupvars = c("Timepoint","Genotype"))
plot <-
  ggplot(datasum, aes(x=Timepoint, y=ME11, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=ME11-se,ymax=ME11+se), width=0.1)+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF")) +
  scale_x_discrete(limits=c("0","p-100","p-600","preBloom-early","p-1000","preBloom-late")) +
  theme_classic(base_size = 24)
plot

ApricotAll$Genotype <- factor(ApricotAll$Genotype, levels = c("A2137", "A1956", "A660","A1267"))
datasum <- summarySE(ApricotAll, measurevar = "ME11", groupvars = c("Timepoint","Genotype"))
plot <- ggplot(datasum, aes(x=Timepoint, y=ME11, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=ME11-se,ymax=ME11+se), width=0.1) + #scale_y_continuous(limits = c(0,1000)) +
  scale_x_discrete(limits=c("a-0","a-100","a-400","Bud-800","Sepal","Petal","Flower-800")) + 
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF"))+
  theme_classic(base_size = 24)
  #theme(text = element_text(size=24),panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot

datAll$Genotype <- factor(datAll$Genotype, levels = c("A209", "A340", "A318", "A323","A2137","A1956","A660","A1267"))
datAll$Timepoint <- gsub("preChill","p-0", datAll$Timepoint)
melt_datAll <- melt(datAll)
plot <-ggplot(melt_datAll, aes(x=Timepoint, y=value, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_smooth(se = FALSE, method = "auto")+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF","#FF9999","#FF3333","#99CCFF","#0033FF")) +
  scale_x_discrete(limits=c("p-0","p-100","p-600","preBloom-early","p-1000","preBloom-late","a-0","a-100","a-400","Bud-800","Sepal","Petal","Flower-800")) +
  theme_classic(base_size = 20)+ facet_wrap( ~ variable, nrow = 4, ncol = 4)+ylab("Eigengene expression")+
  theme(axis.text.x = element_text(angle = 90,size = 12))
ggsave("coMEgenes.png", height = 12, width = 18)

# Select modules
modules = "15";
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleLabels, modules));
modProbes = as.data.frame(probes[inModule])
modProbes$locusName = gsub("\\.v2.1",'',modProbes$`probes[inModule]`)
modProbes_anno =  join(modProbes,annot,by="locusName",type = 'left', match = 'first')
write.csv(modProbes_anno,'co_ME15genes.csv')

for (i in 0:20) {
  # Select modules
  modules = i;
  # Select module probes
  probes = names(datExpr0)
  inModule = is.finite(match(moduleLabels, modules));
  modProbes = as.data.frame(probes[inModule])
  #print(modules)
  print(grep("pp",modProbes$`probes[inModule]`))
}
