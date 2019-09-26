library(reshape2)
library(circlize)

#---------------------circal plot-----------------------------
cytoband.data <- read.table("Pp_chrominfo.txt",header = F, sep = "\t") # peach chromosome range file, convert from gff
cytoband.df <- cytoband.data[,c(1,4,5)]
cytoband.df$V1 <- as.character(cytoband.df$V1)
cytoband.df$V4 <- as.numeric(cytoband.df$V4)
cytoband.df$V5 <- as.numeric(cytoband.df$V5)
circos.initializeWithIdeogram(cytoband.df, plotType = c("axis", "labels"))

Apr_DEgene <- read.csv("Sig_genes_time_1v2_new.csv", header = TRUE, row.names = 1)
Pea_DEgene <- read.csv("../Peach/Sig_genes_600.csv", header = TRUE, row.names = 1)
genelocation <- read.table("geneloc.txt", header = F)
genelocation <- genelocation[,c(1,4,5,9)]

genelocation$V9 <- colsplit(genelocation$V9,";",c('a', 'b'))
genelocation$V9$a <- gsub("ID=","",genelocation$V9$a)
genelocation$V9 <- genelocation$V9[,c(1)]
rownames(genelocation) <- genelocation$V9

apr <- genelocation[ rownames(genelocation) %in% rownames(Apr_DEgene), c(1,2,3)]
peach <- genelocation[ rownames(genelocation) %in% rownames(Pea_DEgene), c(1,2,3)]
apr$V6 <- -(Apr_DEgene$log2FoldChange)
peach$V6 <- Pea_DEgene$log2FoldChange
bed_list = list(apr, peach)

# QTL region:
QTL <- rbind(genelocation[grep("Prupe.1G005500", genelocation$V9):grep("Prupe.1G084200", genelocation$V9),],
            genelocation[grep("Prupe.1G445500", genelocation$V9):grep("Prupe.1G472700", genelocation$V9),],
            genelocation[grep("Prupe.1G529800", genelocation$V9):grep("Prupe.1G534200", genelocation$V9),],
            genelocation[grep("Prupe.2G106200", genelocation$V9):grep("Prupe.2G110900", genelocation$V9),],
            genelocation[grep("Prupe.4G013800", genelocation$V9):grep("Prupe.4G050500", genelocation$V9),],
            genelocation[grep("Prupe.4G152600", genelocation$V9):grep("Prupe.4G205100", genelocation$V9),],
            genelocation[grep("Prupe.5G061500", genelocation$V9):grep("Prupe.5G099600", genelocation$V9),],
            genelocation[grep("Prupe.6G288000", genelocation$V9):grep("Prupe.6G319600", genelocation$V9),],
            genelocation[grep("Prupe.7G148700", genelocation$V9):grep("Prupe.7G210300", genelocation$V9),],
            genelocation[grep("Prupe.8G164300", genelocation$V9):grep("Prupe.8G185000", genelocation$V9),])
bed <- QTL[,c(1,2,3)]
            
# drawing circle plot
circos.initializeWithIdeogram(cytoband.df, plotType = c("axis", "labels"))
circos.genomicRainfall(bed_list, pch = 16, cex = 0.5, col = add_transparency(c("#0033FF", "#FF3333"), transparency = 0.8),bg.border = NA, bg.col = "#FF000019")
circos.genomicDensity(bed_list[[1]], col = c("#0033FF"), track.height = 0.1,bg.border = NA)
circos.genomicDensity(bed_list[[2]], col = c("#FF3333"), track.height = 0.1)
circos.clear()

circos.initializeWithIdeogram(cytoband.df, plotType = c("axis", "labels"))
circos.genomicTrackPlotRegion(bed, 
                              stack = TRUE, 
                              panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "black", border = NA)
}, bg.border = NA, bg.col = "white", track.height = 0.05)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = ifelse(value > 0, "#FF3333", "black"))
                    })
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.5, col = add_transparency(c("#0033FF", "#FF3333"), transparency = 0.8), bg.border = NA, bg.col = "#FFCC0019")
circos.genomicDensity(bed_list[[1]], col = c("#009966"), track.height = 0.1, bg.border = NA)
circos.genomicDensity(bed_list[[2]], col = c("#6666CC"), track.height = 0.1, bg.border = NA)
circos.clear()
