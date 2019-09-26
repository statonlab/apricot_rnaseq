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
  theme_classic()+
  theme(legend.position="none", text = element_text(size = 24))

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

