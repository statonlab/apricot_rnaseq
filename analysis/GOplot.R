setwd("~/Desktop/Jiali/UTK/Apricot/")

library(plyr)
library(tidyr)

summary(Time_1v2)
res_Sig1v2 <- Time_1v2[ which(Time_1v2$padj<0.05),]
write.csv(x = res_Sig1v2, file = "Sig_genes_1v2.csv")

res_Sig1v2 <- read.csv("Sig_genes_1v2.csv",header = TRUE)
res_Sig1v2$X = gsub('\\.v2\\.1', '', res_Sig1v2$X)
colnames(res_Sig1v2)[1] <- "locusName"
match1 <- join(res_Sig1v2, annotat, by = "locusName", match = "first")
match2 <- match1[c(-8,-10)] #remove 8-11th columns in match1, containing names.
#locusname <- read.csv("Ppersica_298_v2.1.locus_transcript_name_map.txt", header = TRUE, comment.char = '', sep='\t')
#colnames(locusname)[1] <- "locusName"
match3 <- join(match2, locusname, by = "locusName", match = "first")
colnames(match3)
match3 <- match3[c(-18,-19)]
head(match3) 
#write.csv(match3, file = "DEgene_400_new.csv")


Gene_up <- match3[which(match3$log2FoldChange>0),]
Gene_down <- match3[which(match3$log2FoldChange<0),]
#head(match3$GO)
Gene_up$GO <- as.character(Gene_up$GO)
s <- strsplit(Gene_up$GO, split = ",")
GO_1v2up <- data.frame(locus = rep(Gene_up$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_1v2up, file = "GO_1v2up", sep = " ", quote = F, row.names = F, col.names = F)

Gene_down$GO <- as.character(Gene_down$GO)
s <- strsplit(Gene_down$GO, split = ",")
GO_1v2down <- data.frame(locus = rep(Gene_down$locusName, sapply(s, length)), GO= unlist(s))
write.table(GO_1v2down, file = "GO_1v2down", sep = " ", quote = F, row.names = F, col.names = F)


############GOplot
library(GOplot)

#---------------AgriGO enrichment---------------------------
sigGenes <- read.table("AgriGo.txt", sep = "\t", header = T)
genes <- sigGenes[,c("GO_acc","entries")]
genes$entries <- gsub("//",",", genes$entries)
s <- strsplit(genes$entries, split = ",")
genes <- data.frame(ID = rep(genes$GO_acc, sapply(s, length)), genes = unlist(s))
genes <- genes[!(genes$genes == ""), ]

sigGenes <- sigGenes[,c(1,2,3,4,9)]
colnames(sigGenes) <- c("ID","category", "term", "count","adj_pval")
GOanalysis <- merge(sigGenes,genes,by="ID", all=T)
GOanalysis$genes <- gsub("\\s", "", GOanalysis$genes)

DEgenes <- read.csv("DEgene_time1v2.csv", header = T)
DEgenes <- DEgenes[,c(1,3)]
names(DEgenes) <- c("genes","logFC")
GOanalysis <- join(GOanalysis,DEgenes, by = "genes", match="first")
GOanalysis$logFC = -GOanalysis$logFC


for (i in 1:9024) {
ID = GOanalysis$ID[i]
GOterm <- GOanalysis[which(GOanalysis$ID == ID),]
up <- sum(GOterm$logFC > 0)
down <- sum(GOterm$logFC < 0)
total <- length(GOterm$logFC)
zscore <- (up-down)/total
GOanalysis$zscore[i] <- zscore
}

GOanalysis <- GOanalysis[c("category","ID","term","count","genes","logFC","adj_pval","zscore")]
GOanalysis$category <- gsub("F","MF",GOanalysis$category)
GOanalysis$category <- gsub("P","BP",GOanalysis$category)
GOanalysis$category <- gsub("C","CC",GOanalysis$category)

GOBubble(GOanalysis, labels = 3)
GOBubble(GOanalysis, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 2, ID=F)  
reduced_GOanalysis <- reduce_overlap(GOanalysis, overlap = 1.6)
GOBubble(reduced_GOanalysis, labels = 3)
GOcir <- GOanalysis[order(GOanalysis$adj_pval),]
sizeGrWindow(15,6)
#pdf("DE_GO_top12.pdf")
GOCircle(GOcir, nsub = 12)
dev.off()
