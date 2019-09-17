library(DESeq2)
library(vegan)
library(ggplot2)

# ANOVA using distance matrices by adonis2 in "vegan"
# calculating distances by rlog
dist_rld <- assay(rld)
# convert data formats
dist_rld_t <- as.data.frame(t(dist_rld))
#head(dist_rld_t)
# create data discription
dist_rld_t.env <- MydataDesign
dist_rld_t.env$Sample <- rownames(MydataDesign)
dist_rld_t.env <- merge(dist_rld_t.env, metrics, by = "Sample")
dist_rld_t.env$rRNA <- dist_rld_t.env$rRNA/100
rownames(dist_rld_t.env) <- dist_rld_t.env$Sample
dist_rld_t.env <- dist_rld_t.env[,-1]
head(dist_rld_t.env)
# adonis2 analysis
adonis2(dist_rld_t ~ Timepoint*Genotype*rRNA, data = dist_rld_t.env)

# contamination removal correction
head(mydata)
head(dist_rld_t.env)
mydata_new <- mydata
for (i in 1:69){
  mydata_new[,i] <- mydata_new[,i] / (1- dist_rld_t.env$rRNA[i])
}
mydata_new <- round(mydata_new, digits = 0)
head(mydata_new)
head(mydata)

# run DESeq again with corrected gene counts
dds_new = DESeqDataSetFromMatrix(countData = mydata_new, colData = MydataDesign, design = ~ phenotype_time)
dds_new <- DESeq(dds_new)
#dds_new = dds_new[rowSums(counts(dds_new))>1, ]
#dds_new <- DESeq(dds_new)
res_new = results(dds_new, contrast = c('phenotype_time', 'LowChill_600', 'HighChill_600'), lfcThreshold = 1, alpha = 0.05)

# plot Padj and corrected Padj
gene_padj <- as.data.frame(res[ which(res$padj <= 1),])[,-c(1,2,3,4)]
gene_padj_new <- as.data.frame(res_new[ which(res_new$padj <= 1),])[,c(5,6)]

Pvalue <- merge(gene_padj, gene_padj_new, by = "row.names", all=T)
rownames(Pvalue) <- Pvalue$Row.names
Pvalue <- Pvalue[,-1]
logPvalue <- -log10(Pvalue)
ggplot(logPvalue)+
  geom_point(mapping = aes(padj.x, padj.y), alpha = 0.1)


# Remove contatmination based on rRNA percentage
data9_env <- dist_rld_t.env[-c(which(dist_rld_t.env$rRNA > 0.80)),]
data9 <- dist_rld_t[-c(which(dist_rld_t.env$rRNA > 0.50)),]
adonis2(data9 ~ Timepoint*Genotype*rRNA, data = data9_env)
