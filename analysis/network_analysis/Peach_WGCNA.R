setwd("~/Desktop/Jiali/UTK/Peach//")
## WCNNA co-expression analysis
source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#install.packages("WGCNA")
library(DESeq2)
library(WGCNA)
options(stringsAsFactors = F)
PeachData = read.csv("normalized Peach.csv",row.names = 1, header = T)
dim(PeachData)
names(PeachData)
datExpr0 = as.data.frame(t(PeachData))
datExpr0 <- as.data.frame(lapply(datExpr0,as.numeric)) #datExpr0 has to be numeric.
rownames(datExpr0) <- colnames(PeachData)
datExpr0[1:69,1:3]
# column: gene names, row: sample variables

# check missing value and outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3)
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

# cluster the samples
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


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
net = blockwiseModules(datExpr0, power = 8,
                       TOMType = "unsigned", minModuleSize = 300,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PeachTOM",
                       verbose = 3)
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs1 = net$MEs;
geneTree = net$dendrograms[[1]];

# Load genotype data
datTraits0 = read.csv("PeachTrait.csv", header = T,row.names = 1,stringsAsFactors = T)
rownames(datTraits0) 
rownames(datExpr0)

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits0[,c(10,11,12)], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(11,9)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#pdf("corrplot_0116.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0)[10:12],
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

#---------------------------------------------------

#==========visualize gene network============================================
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 8);
# Read in the annotation file
annot = read.csv(file = "../Apricot/Ppersica_298_v2.1.annotation_info.txt",header = T,sep = "\t");
annot$arabi.hit = gsub('.{2}$',"",annot$Best.hit.arabi.name)
annot$arabi.hit = gsub("T","t",annot$arabi.hit)
annot$arabi.hit = gsub("G","g",annot$arabi.hit)

## export to cytoscape
# Select modules
modules = "5";
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
                               edgeFile = paste("new_CytoscapeInput-edges-ME2.txt", sep=""),
                               nodeFile = paste("new_CytoscapeInput-nodes-ME2.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# DEgenes from DEseq
DEprof = read.csv("DEgene_600.csv")
DEgenes = DEprof$locusName
# Genes in each module
ModuleTab <- as.data.frame(matrix(NA,nrow=21,ncol=2))
names(ModuleTab) <- c("Non-DEgenes","DEgenes")
row.names(ModuleTab) <- names(MEs1)
for (i in 0:20) {
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

pdf( "MEbar_0327.pdf", width = 7, height = 6 )
par(mfrow = c(1,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0)[10:12],
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
        horiz = T, legend = names(ModuleTab),args.legend = list(x="right"))
dev.off()

# Expression profile of each module
# 1. the genes inside a given module are summarize with the module eigengene, which can be considered as the best summary of the standardized module expression data
datAll <- cbind(MEs1,MydataDesign[,c(1,3)])
datAll$Genotype <- factor(datAll$Genotype, levels = c("A209", "A340", "A318", "A323"))
datAll$Timepoint <- gsub("preChill","0", datAll$Timepoint)
datasum <- summarySE(datAll, measurevar = "ME4", groupvars = c("Timepoint","Genotype"))
plot <-
  ggplot(datasum, aes(x=Timepoint, y=ME4, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=ME4-se,ymax=ME4+se), width=0.1)+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF")) +
  scale_x_discrete(limits=c("0","100","600","preBloom-early","1000","preBloom-late")) +
  theme_classic(base_size = 24)
plot

melt_datAll <- melt(datAll)
plot <-ggplot(melt_datAll, aes(x=Timepoint, y=value, color=Genotype, group=Genotype)) +
  geom_point(size=2) + geom_smooth(se = FALSE, method = "auto")+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF")) +
  scale_x_discrete(limits=c("0","100","600","preBloom-early","1000","preBloom-late")) +
  theme_classic(base_size = 20)+ facet_wrap( ~ variable, nrow = 3, ncol = 7)+ylab("Eigengene expression")+
  theme(axis.text.x = element_text(angle = 90,size = 12))
ggsave("peachMEgenes.png", height = 10, width = 20)

#---------------------20190304-----------------------
# Select modules
modules = "5";
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleLabels, modules));
modProbes = as.data.frame(probes[inModule])
modProbes$locusName = gsub("\\.v2.1",'',modProbes$`probes[inModule]`)
modProbes_anno =  join(modProbes,annotat,by="locusName",type = 'left', match = 'first')
write.csv(modProbes_anno,'new_ME5genes.csv')

OverlapList <- as.data.frame(matrix(NA,nrow=21,ncol=3))
names(OverlapList) <- c('ME','X','T')
for (i in 0:20) {
  modules = as.character(i)
  inModule = is.finite(match(moduleLabels, modules))
  modProbes = gsub("\\.v2.1",'',probes[inModule])
  OverlapList[i+1,1] = modules
  OverlapList[i+1,2] = length(intersect(Apri_ME15$locusName,modProbes))
  OverlapList[i+1,3] = length(modProbes)
}
plot.new()
plot(OverlapList$ME,OverlapList$X,type = "o",col = "red", xlab = "MEs",ylab = "Number of genes", 
     main = "Overlapped genes with Apricot ME15")
lines(OverlapList$X, type = "o", col = "blue")


