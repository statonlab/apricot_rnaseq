setwd("~/Desktop/Jiali/UTK/Apricot/")
## WCNNA co-expression analysis
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library(WGCNA)
options(stringsAsFactors = F)
ApricotData = read.csv("normalizedData.csv")
#dim(ApricotData)
names(ApricotData)
datExpr0 = as.data.frame(t(ApricotData))
names(datExpr0) = ApricotData$X
datExpr0 <- as.data.frame(lapply(datExpr0,as.numeric)) #datExpr0 has to be numeric.
rownames(datExpr0) = names(ApricotData)
datExpr0 = datExpr0[-1,]
rownames(datExpr0)= gsub("X","",rownames(datExpr0))
datExpr0[1:3,1:3]
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

# Plot a line to show the cut
abline(h = 400000, col = "red"); # I don't want to lose samples, so I skip the following steps
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Auto one step construction
enableWGCNAThreads() # Not work in Rstudio, this produce an error in the softThread step, need to change to allowWGCNAThreads()
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
net = blockwiseModules(datExpr0, power = 11,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ApricotTOM",
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
# save(MEs, moduleLabels, moduleColors, geneTree,
     #file = "FemaleLiver-02-networkConstruction-auto.RData")

# Load genotype data
datTraits0 = read.csv("ApricotTraits.csv", header = T,row.names = 1,stringsAsFactors = T)
rownames(datTraits0) = rownames(datExpr0)
datTraits1 = datTraits0[,5:9]

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits1, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#pdf("corrplot.pdf")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0),
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

#----------------------------------------------------
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



#==========visualize gene network============================================
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 11);
# Read in the annotation file
annot = read.csv(file = "Ppersica_298_v2.1.annotation_info.txt",header = T,sep = "\t");
annot$arabi.hit = gsub('.{2}$',"",annot$Best.hit.arabi.name)
annot$arabi.hit = gsub("T","t",annot$arabi.hit)
annot$arabi.hit = gsub("G","g",annot$arabi.hit)
# Select module probes, extract genes in one module
probes = names(datExpr0)
inModule = (moduleColors=="midnightblue");
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )

nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )

## export to cytoscape
# Select modules
modules = "20";
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
                               edgeFile = paste("CytoscapeInput-edges-ME20.txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-ME20.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


# DEgenes from DEseq
DEprof = read.csv("DEgene_time1v2.csv")
DEgenes = DEprof$locusName
# Genes in each module
ModuleTab <- as.data.frame(matrix(NA,nrow=23,ncol=2))
names(ModuleTab) <- c("Non-DEgenes","DEgenes")
row.names(ModuleTab) <- names(MEs1)
for (i in 0:22) {
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
pdf("MEbar.pdf")

sizeGrWindow(11,8)
par(mfrow = c(1,2))

barplot(dat, main="Number of genes in each module",
        xlab="Number of genes", las=2, col=c("blue","red"),
        horiz = T, legend = names(ModuleTab))
dev.off()

# Expression profile of each module
# 1. the genes inside a given module are summarize with the module eigengene, which can be considered as the best summary of the standardized module expression data
ExptDesign <- read.csv("MydataDesign_1.csv",header = T,row.names = 1)
datAll <- cbind(MEs1,ExptDesign[,c(1,3)])
datAll$Stage <- gsub("C","Sepal", datAll$Stage)
datAll$Stage <- gsub("D","Petal", datAll$Stage)
datAll$Genotype <- factor(datAll$Genotype, levels = c("A2137","A1956","A660", "A1267"))
datasum <- summarySE(datAll, measurevar = "ME2", groupvars = c("Stage","Genotype"))
plot16 <- 
  ggplot(datasum, aes(x=Stage, y=ME2, color=Genotype, group=Genotype))+
  geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=ME2-se,ymax=ME2+se), width=0.1)+ 
  theme(text = element_text(size=6), axis.text = element_text(size = 6)) + 
  theme_classic(base_size = 20)+scale_x_discrete(limits=c("0","100","400","Bud-800","Sepal","Petal","Flower-800"))+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF"))
plot16

melt_datAll <- melt(datAll)
plot <-ggplot(melt_datAll, aes(x=Stage, y=value, color=genotype, group=genotype)) +
  geom_point(size=2) + geom_smooth(se = FALSE, method = "auto")+
  scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF")) +
  scale_x_discrete(limits=c("0","100","400","Bud-800","Sepal","Petal","Flower-800")) +
  theme_classic(base_size = 20)+ facet_wrap( ~ variable, nrow = 4, ncol = 6)+ylab("Eigengene expression")+
  theme(axis.text.x = element_text(angle = 90,size = 12))
ggsave("apricotMEgenes.png", height = 12, width = 16)

# Select modules
library(plyr)
modules = "20";
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleLabels, modules));
modProbes = as.data.frame(probes[inModule])
modProbes$locusName = gsub("\\.v2.1",'',modProbes$`probes[inModule]`)
modProbes_annot = join(modProbes,annot,by = 'locusName',match='first')
write.csv(modProbes_annot,'ME20_genes.csv')

for (i in 0:22) {
  modules = i;
  # Select module probes
  probes = names(datExpr0)
  inModule = is.finite(match(moduleLabels, modules));
  modProbes = as.data.frame(probes[inModule])
  print(grep('pp', modProbes$`probes[inModule]`))
}


