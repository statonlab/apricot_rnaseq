setwd("~/Desktop/Jiali/UTK/Peach/")
filenames <- list.files("conetwork", pattern = "*.txt", full.names = T)
ldf <- lapply(filenames, function(x) read.table(x, header = T, sep = "\t"))
GOterms <- do.call(rbind, ldf)[,1:3]
GOterms <- GOterms[!duplicated(GOterms$GO_acc),]
BP <- GOterms[which(GOterms$term_type == "P"),]
MF <- GOterms[which(GOterms$term_type == "F"),]
CC <- GOterms[which(GOterms$term_type == "C"),]
MEs <- gsub("conetwork/","",filenames)
MEs <- gsub(".txt","",MEs)
for (i in 1:16) {
  BP <- merge(BP, as.data.frame(ldf[i])[1:5,c(1,9)], all.x=TRUE)
  names(BP)[i+3] <- MEs[i]
}

BP_Total <- as.data.frame(matrix(NA,nrow=16,ncol=1))
rownames(BP_Total) <- MEs
for (i in 1:16) {
  c <- length(which(BP[,i+3]<0.05))
  BP_Total[i,1] <-c
}

MEs1 <- c("ME2","ME10","ME14","ME8","ME15","ME1","ME7","ME3","ME6","ME9","ME5","ME12","ME13","ME4","ME11","ME0")
BP[is.na(BP)] <- 1
filterBP <- BP[rowSums(BP[,-c(1:3)] < 0.05) >= 1, ]
melt_BP_filtered <- melt(data = filterBP, id.vars = "Term", measure.vars = MEs)

p1 <- ggplot(melt_BP_filtered,aes(x=variable,y=Term,fill=value))+
  geom_tile()+
  labs(fill='FDR') +
  scale_x_discrete(limits=rev(MEs1))+
  scale_fill_gradient(low = "red", high = "white", limits = c(0,1)) +
  labs(x = "",y = "GO Term") +coord_flip()+
  theme(axis.text.x = element_text(angle=75,hjust=1,vjust=1.0, size = 16),
        axis.text.y = element_text(size = 14))

rownames(BP_Total) <- factor(rownames(BP_Total), levels = MEs1)
p <- ggplot(data=BP_Total, aes(x=rownames(BP_Total), y=V1)) +
  geom_bar(stat="identity", fill="red") + coord_flip() + scale_x_discrete(limits=rev(MEs1))+
  theme(text = element_text(size=14),panel.background = element_blank(), axis.text.y = element_text(size = 12,colour = "black"),axis.line.y = element_blank(), axis.ticks.y = element_blank())
p2<- p + labs(title = "Number of significant GO terms", x = "", y = "Number of terms") +geom_text(aes(label = paste0("n=",V1), y = V1 + 5),
                                                            position = position_dodge(0.9),
                                                            vjust = 0.4,
                                                            color   = "black")

p2
grid.newpage()
ggsave("conetwork/GObar.png", width = 18, height = 12)
