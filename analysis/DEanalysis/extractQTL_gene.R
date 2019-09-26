setwd("~/Desktop/Jiali/UTK/Apricot/")
library(stringi)
library("VennDiagram")
library(DESeq2)

# Overlap genes between DE genes and QTL genes:
transcripts <- rownames(mydata)
DEgene1v2 <- rownames(resTime_1v2) # DE genes between endo and eco

DEGene_Endo <- c(rownames(read.csv("Sig_genes_0.csv", header = TRUE, row.names = 1)),
                 rownames(read.csv("Sig_genes_100.csv", header = TRUE, row.names = 1)),
                 rownames(read.csv("Sig_genes_400.csv", header = TRUE, row.names = 1)))
DEgene_inter <- rownames(read.csv("../Peach/Sig_genes_inter.csv", header = T, row.names = 1)) 

DEgene_peach <- rownames(read.csv("../Peach/Sig_genes_600.csv", header = T, row.names = 1))

find_inter <- function(gene1,gene2) {
  QTLgene <- transcripts[grep(gene1, transcripts):grep(gene2, transcripts)]
  print(length(QTLgene))
  intersect(DEgene_peach,QTLgene)
  
  #print(length(intersect(DEgene1v2,QTLgene)))
}

# region qCR1d-2008
find_inter("Prupe.1G005500","Prupe.1G084200")

# region qCR1c-2009
find_inter("Prupe.1G445500","Prupe.1G472700")

# region qCR1a-2009
find_inter("Prupe.1G529800","Prupe.1G534200")

# region qCR2-2009
find_inter("Prupe.2G106200","Prupe.2G110900")

# region qCR4a-2008/2009
find_inter("Prupe.4G013800", "Prupe.4G050500")

#region qCR4b-2008/2009
find_inter("Prupe.4G152600","Prupe.4G205100")

# region qCR5-2008/2009
find_inter("Prupe.5G061500","Prupe.5G099600")

# region qCR6-2008
find_inter("Prupe.6G288000","Prupe.6G319600")

# region qCR7-2008/2009
find_inter("Prupe.7G148700","Prupe.7G210300")

# region qCR8-2008/2009
find_inter("Prupe.8G164300","Prupe.8G185000")

# bud break
# region qBD1c-2007 and qBD1b-2006/2007
find_inter("Prupe.1G161000", "Prupe.1G379300") # couldn't find one marker, a large region

# region qBD2-2009

# region qBD3-2008
find_inter("Prupe.3G045100","Prupe.3G051100")

# region qBD5-2008

# region qBD7b-2007
find_inter("Prupe.7G067900","Prupe.7G082900")

# region qBD8-2008
find_inter("Prupe.8G029800", "Prupe.8G036600")
