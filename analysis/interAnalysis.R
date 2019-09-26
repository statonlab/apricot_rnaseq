setwd("~/Desktop/Jiali/UTK/Peach/")
library(ggplot2)

## Compare the DEGs of apricot and peach
apricot_DE <- read.csv('../Apricot/Sig_genes_time_1v2_new.csv',header = TRUE)
peach_DE <- read.csv('Sig_genes_600.csv', header = TRUE)
apricot_DE$log2FoldChange = -apricot_DE$log2FoldChange
all <- as.data.frame(unique(c(as.vector(apricot_DE$X),as.vector(peach_DE$X))))
names(all) <- 'X'
apricot_results <- read.csv('apricot_results.csv',header = T)
apricot_results$log2FoldChange = -apricot_results$log2FoldChange
peach_results <- read.csv('DE_results.csv',header = T)
all <- join(all, apricot_results[,c(1,3,7)],by = 'X',type='left')
names(all) <- c('X','ApricotFC','ApricotPvalue')
all <- join(all,peach_results[,c(1,3,7)], by='X',type='left')
all$ApricotPvalue[which(is.na(all$ApricotPvalue))] = 1
all$padj[which(is.na(all$padj))] = 1
all$Significant <- 'Both'
all$Significant[which(all$ApricotPvalue<0.05 & all$padj>0.05)] <- 'Apricot'
all$Significant[which(all$ApricotPvalue>0.05 & all$padj<0.05)] <- 'Peach'
all_annot <- all
all_annot$X = gsub('\\.v2\\.1', '', all_annot$X)
colnames(all_annot)[1] <- "locusName"
annot <- read.csv("../Apricot/Ppersica_298_v2.1.annotation_info.txt", header = TRUE, comment.char = '', sep='\t')
match1 <- join(all_annot, annot, by = 'locusName', match = "first")
all_annot <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
write.csv(all_annot,'DE_all.csv')

#--------------Volcano plot with DEgenes in peach and apricot-----------------------------
apricot_results$logFDR <- -log10(apricot_results$padj)
apricotData <- apricot_results[-which(apricot_results$logFDR == 0),]
apricotData$threshold = as.factor(apricotData$padj <= 0.05 & apricotData$log2FoldChange >=1 | apricotData$padj<= 0.05 & apricotData$log2FoldChange <= (-1))


ggplot(data = apricotData, aes(x=log2FoldChange, y=logFDR,colour=threshold)) +
  geom_point(size=1, alpha = 0.5) + scale_color_manual(values = c("black","red"))+
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 85),colour='black',linetype=2)+
  geom_segment(aes(x = -15, y = 1.3, xend = 10, yend = 1.3),colour='black',linetype=2)+
  geom_segment(aes(x = -1, y = 0, xend = -1, yend = 85),colour='black',linetype=2)+
  labs( x = "log2 fold change", y="-log10 (adjusted p-value)")+
  theme_classic(base_size = 20)+
  theme(legend.position="none")
ggsave("volcano.png", height = 5, width = 5)

peach_results$logFDR <- -log10(peach_results$padj)
peachData <- peach_results[-which(peach_results$logFDR == 0),]
peachData$threshold <- as.factor(peachData$padj <=0.05 & peachData$log2FoldChange >=1 | peachData$padj <=0.05 & peachData$log2FoldChange <= -1)
ggplot(data = peachData, aes(x=log2FoldChange, y=logFDR,colour=threshold)) +
  geom_point(size=1, alpha = 0.5) + scale_color_manual(values = c("black","red"))+
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 85),colour='black',linetype=2)+
  geom_segment(aes(x = -10, y = 1.3, xend = 10, yend = 1.3),colour='black',linetype=2)+
  geom_segment(aes(x = -1, y = 0, xend = -1, yend = 85),colour='black',linetype=2)+
  labs( x = "log2 fold change", y="-log10 (adjusted p-value)")+
  theme_classic()+
  theme(legend.position="none", text = element_text(size = 24))

#Generate files for GO enrichment
Apri_only <- all_annot[which(all_annot$Significant == 'Apricot'),]
Peach_only <- all_annot[which(all_annot$Significant=='Peach'),]
Apri_only_At <- data.frame(gsub('\\.\\d','', Apri_only$Best.hit.arabi.name))
Peach_only_At <- data.frame(gsub('\\.\\d','',Peach_only$Best.hit.arabi.name))
write.table(Apri_only_At,'ApricotOnly_AT.txt',sep = "\n", quote = F, row.names = F, col.names = F)
write.table(Peach_only_At,'PeachOnly_AT.txt',sep = "\n", quote = F, row.names = F, col.names = F)

Apri_only$GO <- as.character(Apri_only$GO)
s <- strsplit(Apri_only$GO, split = ",")
GO_AprOnly <- data.frame(locus = rep(Apri_only$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_AprOnly, file = "GO_Apri_only", sep = " ", quote = F, row.names = F, col.names = F)

Peach_only$GO <- as.character(Peach_only$GO)
s <- strsplit(Peach_only$GO, split = ",")
GO_PeachOnly <- data.frame(locus = rep(Peach_only$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_PeachOnly, file = "GO_Peach_only", sep = " ", quote = F, row.names = F, col.names = F)

#-----Plot all DE genes------------------------------------
ggplot(all, aes(x=ApricotFC, y=log2FoldChange, color=Significant)) +
  geom_point(size=1, alpha=0.5)+
  geom_segment(aes(x = 1, y = 1, xend = 1, yend = 10),colour='black',linetype=2)+
  geom_segment(aes(x = 1, y = 1, xend = 10, yend = 1),colour='black',linetype=2)+
  geom_segment(aes(x = -1, y = -1, xend = -1, yend = -10),colour='black',linetype=2)+
  geom_segment(aes(x = -1, y = -1, xend = -10, yend = -1),color='black',linetype=2)+
  xlab('logFC in Apicot') + ylab('logFC in Peach')+
  scale_color_brewer(palette="Dark2")+
  theme_classic(base_size = 16)
  #scale_color_manual(values=c('#0033FF', '#00CC99','#FF3333'))

## Compare the endo and eco clusters from apricot and peach
Apri_ME2 <- read.csv('../Apricot/ME2_genes.csv',header = TRUE,row.names = 1)
Apri_ME15 <- read.csv('../Apricot/ME15_genes.csv',header = TRUE,row.names = 1)
Apri_ME20 <- read.csv('../Apricot/ME20_genes.csv', header = TRUE, row.names = 1)
#Peach_ME15 <- read.csv('ME15genes.csv', header = T,row.names = 1)
#Peach_ME19 <- read.csv('ME19genes.csv', header = T,row.names = 1)
#Peach_ME20 <- read.csv('ME20genes.csv',header = T,row.names = 1)
#Peach_ME6 <- read.csv('ME6genes.csv',header = T,row.names = 1)
Peach_6 <- read.csv('new_ME6genes.csv',header = T,row.names = 1)
Peach_4 <- read.csv('new_ME4genes.csv',header = T,row.names = 1)
Peach_10 <- read.csv('new_ME10genes.csv',header = T,row.names = 1)
Peach_14 <- read.csv('new_ME14genes.csv',header = T,row.names = 1)
Peach_15 <- read.csv('new_ME15genes.csv',header = T,row.names = 1)
#Eco
length(intersect(Apri_ME15$locusName,Peach_4$locusName))
length(intersect(Apri_ME15$locusName,Peach_10$locusName))
length(intersect(Apri_ME15$locusName,Peach_14$locusName))
length(intersect(Apri_ME15$locusName,Peach_15$locusName))
length(intersect(Apri_ME15$locusName,Peach_6$locusName))

length(intersect(Apri_ME20$locusName,Peach_4$locusName))

#Endo

length(intersect(Apri_ME2$locusName,Peach_6$locusName))
length(intersect(Apri_ME2$locusName,Peach_14$locusName))
length(intersect(Apri_ME2$locusName,Peach_15$locusName))
length(intersect(Apri_ME2$locusName,Peach_4$locusName))
length(intersect(Apri_ME2$locusName,Peach_10$locusName))

length(Apri_ME15$locusName)
length(Apri_ME2$locusName)
length(Peach_6$locusName)
length(Peach_4$locusName)
length(Peach_10$locusName)
length(Peach_14$locusName)
length(Peach_15$locusName)

# Venn diagram
## ---ones
area2=length(Apri_ME2$locusName)
area1=length(Peach_6$locusName)
area3=length(Peach_10$locusName)

#---pairs
n12=length(intersect(Peach_6$locusName, Apri_ME2$locusName))
n13=length(intersect(Peach_4$locusName,Peach_10$locusName))
n23=length(intersect(Apri_ME15$locusName,Peach_10$locusName))

#---trios
n123=length(Reduce(intersect,list(Apri_ME15,Peach_4,Peach_10)))

grid::grid.newpage()
draw.pairwise.venn(area1, area2, n12,
                 scaled=TRUE, 
                 fill=c("red", "yellow"),
                 cat.cex = 2,
                 cex = 2)
#----------------------------------20190306-----------------------------------
#GO enrichment
install.packages('GOplot')
library(GOplot)

sigGenes <- read.table("GO_inter_enrichment.txt", sep = "\t", header = T)
genes <- sigGenes[,c("GO_acc","entries")]
genes$entries <- gsub("//",",", genes$entries)
s <- strsplit(genes$entries, split = ",")
genes <- data.frame(ID = rep(genes$GO_acc, sapply(s, length)), genes = unlist(s))
genes <- genes[!(genes$genes == ""), ]

sigGenes <- sigGenes[,c(1,2,3,4,9)]
colnames(sigGenes) <- c("ID","category", "term", "count","adj_pval")
GOanalysis <- merge(sigGenes,genes,by="ID", all=T)
GOanalysis$genes <- gsub("\\s", "", GOanalysis$genes)

DEgenes <- read.csv("DEgene_inter.csv", header = T,row.names = 1)
DEgenes <- DEgenes[,c(1,3)]
names(DEgenes) <- c("genes","logFC")
GOanalysis <- join(GOanalysis,DEgenes, by = "genes", match="first")
#GOanalysis$logFC = -GOanalysis$logFC


for (i in 1:4579) {
  ID = GOanalysis$ID[i]
  GOterm <- GOanalysis[which(GOanalysis$ID == ID),]
  up <- sum(GOterm$logFC > 0)
  down <- sum(GOterm$logFC < 0)
  total <- length(GOterm$logFC)
  zscore <- (up-down)/sqrt(total)
  GOanalysis$zscore[i] <- zscore
}

GOanalysis <- GOanalysis[c("category","ID","term","count","genes","logFC","adj_pval","zscore")]
GOanalysis$category <- gsub("F","MF",GOanalysis$category)
GOanalysis$category <- gsub("P","BP",GOanalysis$category)
GOanalysis$category <- gsub("C","CC",GOanalysis$category)
Enriched_GO <- as.vector(unique(GOanalysis$term[which(GOanalysis$adj_pval<0.05)]))
#EC_gene <- GOanalysis[which(GOanalysis$adj_pval<0.05), c('genes','logFC')]
chord <- chord_dat(GOanalysis, DEgenes, Enriched_GO)

GOBubble(GOanalysis, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 2, ID=F)  
reduced_GOanalysis <- reduce_overlap(GOanalysis, overlap = 1.6)
GOBubble(reduced_GOanalysis,labels = 1.3)


GOCluster(GOanalysis, Enriched_GO, clust.by = 'term',lfc.min = -10,lfc.max = 10)
#ggsave('GO_inter.png', height = 7, width = 14)

intersetTab <- as.data.frame(matrix(NA,nrow=21,ncol=23))
names(intersetTab) <- names(MEs1)
row.names(intersetTab) <- names(MEsP)
probes = names(datExpr0)
probesP = names(datExprP)
moduleLabelsP = net_peach$colors
for (i in 0:22) {
  inModule = (moduleLabels==i)
  modProbes = probes[inModule]
  MEc= which(substring(names(MEs1),3)==i)
  for (j in 0:20) {
    inModule_P = (moduleLabelsP==j)
    modProbes_P = probesP[inModule_P]
    MEr = which(substring(names(MEsP),3)==j)
    intersetTab[MEr, MEc] = length(intersect(modProbes, modProbes_P))
  }
}
overlapnumber = as.matrix(intersetTab)
colfunc <- colorRampPalette(c("white", "red"))
sizeGrWindow(11,12)
labeledHeatmap(Matrix = intersetTab,
               xLabels = names(intersetTab),
               yLabels = row.names(intersetTab),
               colorLabels = FALSE,
               colors = colfunc(50),
               textMatrix = overlapnumber,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(0,500),
               main = paste("Number of genes overlapped between apricot and peach"))

#------------------------20190328-------------------GO heatmap with co-expression modules
peachGO_4 <- read.table("AgriGO_ME4.txt", sep = "\t", header = T)
peachGO_10 <- read.table("AgriGO_ME10.txt", sep = "\t", header = T)
peachGO_6 <- read.table("AgriGO_ME6.txt", sep = "\t", header = T)
apricoGO_2 <- read.table("Apricot_ME2.txt", sep = '\t', header = T)
apricoGO_15 <- read.table("Apricot_ME15.txt", sep = '\t', header = T)

all_GO <- rbind(peachGO_10[,1:3], peachGO_4[,1:3], peachGO_6[,1:3], apricoGO_15[,1:3], apricoGO_2[,1:3])
all_GO <- all_GO[!duplicated(all_GO$GO_acc),]
BP <- all_GO[which(all_GO$term_type == "P"),]
MF <- all_GO[which(all_GO$term_type == "F"),]
CC <- all_GO[which(all_GO$term_type == "C"),]

BP <- merge(BP, apricoGO_2[,c(1,9)], all.x=TRUE)
names(BP)[4] <- "Apricot_ME2"
BP <- merge(BP, peachGO_6[,c(1,9)], all.x=TRUE)
names(BP)[5] <- "Peach_ME6" 
BP <- merge(BP, apricoGO_15[,c(1,9)], all.x=TRUE)
names(BP)[6] <- "Apricot_ME15"
BP <- merge(BP, peachGO_4[,c(1,9)], all.x=TRUE)
names(BP)[7] <- "Peach_ME4"
BP <- merge(BP, peachGO_10[,c(1,9)], all.x=TRUE)
names(BP)[8] <- "Peach_ME10"
BP[is.na(BP)] <- 1

BP_filtered <- BP[rowSums(BP < 0.00001) >=1, ]
melt_BP_filtered <- melt(data = BP_filtered, id.vars = "Term", measure.vars = c("Apricot_ME2", "Apricot_ME15","Peach_ME4","Peach_ME6","Peach_ME10"))
ggplot(melt_BP_filtered,aes(x=variable,y=Term,fill=value))+
  geom_tile()+
  labs(fill='FDR') +
  scale_x_discrete(limits=c("Apricot_ME2","Peach_ME6","Apricot_ME15","Peach_ME4","Peach_ME10"))+
  scale_fill_gradient(low = "blue", high = "white", limits = c(0,1)) +
  labs(x = "",y = "GO Term") +
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1.0, size = 12),
        axis.text.y = element_text(size = 10))

MF <- merge(MF, apricoGO_2[,c(1,9)], all.x=TRUE)
names(MF)[4] <- "Apricot_ME2"
MF <- merge(MF, peachGO_6[,c(1,9)], all.x=TRUE)
names(MF)[5] <- "Peach_ME6" 
MF <- merge(MF, apricoGO_15[,c(1,9)], all.x=TRUE)
names(MF)[6] <- "Apricot_ME15"
MF <- merge(MF, peachGO_4[,c(1,9)], all.x=TRUE)
names(MF)[7] <- "Peach_ME4"
MF <- merge(MF, peachGO_10[,c(1,9)], all.x=TRUE)
names(MF)[8] <- "Peach_ME10"
MF[is.na(MF)] <- 1

MF_filtered <- MF[rowSums(MF < 0.00001) >=1, ]
melt_MF_filtered <- melt(data = MF_filtered, id.vars = "Term", measure.vars = c("Apricot_ME2", "Apricot_ME15","Peach_ME4","Peach_ME6","Peach_ME10"))
ggplot(melt_MF_filtered,aes(x=variable,y=Term,fill=value))+
  geom_tile()+
  labs(fill='FDR') +
  scale_x_discrete(limits=c("Apricot_ME2","Peach_ME6","Apricot_ME15","Peach_ME4","Peach_ME10"))+
  scale_fill_gradient(low = "blue", high = "white", limits = c(0,1)) +
  labs(x = "",y = "GO Term") +
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1.0, size = 12),
        axis.text.y = element_text(size = 10))

CC <- merge(CC, apricoGO_2[,c(1,9)], all.x=TRUE)
names(CC)[4] <- "Apricot_ME2"
CC <- merge(CC, peachGO_6[,c(1,9)], all.x=TRUE)
names(CC)[5] <- "Peach_ME6" 
CC <- merge(CC, apricoGO_15[,c(1,9)], all.x=TRUE)
names(CC)[6] <- "Apricot_ME15"
CC <- merge(CC, peachGO_4[,c(1,9)], all.x=TRUE)
names(CC)[7] <- "Peach_ME4"
CC <- merge(CC, peachGO_10[,c(1,9)], all.x=TRUE)
names(CC)[8] <- "Peach_ME10"
CC[is.na(CC)] <- 1

CC_filtered <- CC[rowSums(CC < 0.00001) >=1, ]
melt_CC_filtered <- melt(data = CC_filtered, id.vars = "Term", measure.vars = c("Apricot_ME2", "Apricot_ME15","Peach_ME4","Peach_ME6","Peach_ME10"))
ggplot(melt_CC_filtered,aes(x=variable,y=Term,fill=value))+
  geom_tile()+
  labs(fill='FDR') +
  scale_x_discrete(limits=c("Apricot_ME2","Peach_ME6","Apricot_ME15","Peach_ME4","Peach_ME10"))+
  scale_fill_gradient(low = "blue", high = "white", limits = c(0,1)) +
  labs(x = "",y = "GO Term") +
  theme(axis.text.x = element_text(angle=30,hjust=1,vjust=1.0, size = 12),
        axis.text.y = element_text(size = 10))

#----------------------20190404--------------------------------
sigGenes_Both <- read.csv("BothGO_At.txt", header = T, sep = '\t')
sigGene_Apri <- read.csv("ApriOnlyGO_At.txt", header = T, sep = '\t')
sigGene_Peach <- read.csv("PeachOnly_AT.txt", header = T, sep = '\t')

sigGenes_Both$logFDR <- log10(sigGenes_Both$FDR)

p<-ggplot(data=sigGenes_Both[1:10,], aes(x=reorder(Term, -logFDR), y=-logFDR)) +
    geom_bar(stat="identity", fill="#D95F02", width = 0.8) + coord_flip()+ theme(text = element_text(size=14),panel.background = element_blank(), axis.text.y = element_text(size = 16,colour = "black"),axis.line.y = element_blank(), axis.ticks.y = element_blank())
p + labs(x = "", y = "-log10(adjusted p-value)") +geom_text(aes(label = paste0("n=",queryitem), y = -logFDR + 1),
                                                            position = position_dodge(0.9),
                                                            vjust = 0,
                                                            color   = "black")

sigGene_Apri$logFDR <- log10(sigGene_Apri$FDR)
p<-ggplot(data=sigGene_Apri[1:10,], aes(x=reorder(Term, -logFDR), y=-logFDR)) +
  geom_bar(stat="identity", fill="#1B9E77",width = 0.8) + coord_flip() + theme(text = element_text(size=14),panel.background = element_blank(), axis.text.y = element_text(size = 16,colour = "black"),axis.line.y = element_blank(), axis.ticks.y = element_blank())
p + labs(x = "", y = "-log10(adjusted p-value)") +geom_text(aes(label = paste0("n=",queryitem), y = -logFDR + 0.5),
                                                            position = position_dodge(0.9),
                                                            vjust = 0,
                                                            color   = "black")
sigGene_Peach$logFDR <- log10(sigGene_Peach$FDR)
p<-ggplot(data=sigGene_Peach[1:10,], aes(x=reorder(Term, -logFDR), y=-logFDR)) +
  geom_bar(stat="identity", fill="#7570B3", width = 0.8) + coord_flip() + theme(text = element_text(size=14),panel.background = element_blank(), axis.text.y = element_text(size = 16,colour = "black"),axis.line.y = element_blank(), axis.ticks.y = element_blank())
p + labs(x = "", y = "-log10(adjusted p-value)") +geom_text(aes(label = paste0("n=",queryitem), y = -logFDR + 1.5),
                                                            position = position_dodge(0.9),
                                                            vjust = 0,
                                                            color   = "black")

