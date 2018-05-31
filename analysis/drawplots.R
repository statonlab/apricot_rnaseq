install.packages("stringi")
library(stringi)

geneList <- read.delim("Gene modified list.txt")
stri_sub(geneList$locusName, 15, 2) <- ".v2.1"
geneList <- geneList[1:84,]
head(geneList$locusName)


TimePlot <- function(x) {
  data <- plotCounts(dds, gene=x,intgroup=c("Stage","genotype"), returnData=TRUE)
  plot <-
    ggplot(data, aes(x=Stage, y=count, color=genotype, group=genotype)) +
    geom_point(size=2) + geom_smooth(se = FALSE, method = "loess") +
    scale_y_continuous(limits = c(0.01,10000), trans = "log10") + 
    theme(text = element_text(size=10), 
               plot.margin = unit(c(2, 2, 2, 2), "cm")) + 
    labs(title = x) + scale_color_manual(values=c("#0033FF","#FF3333","#FF9999","#99CCFF"))
  print(plot)
  #  return()
}


pdf('GeneList_plot.pdf')
for (name in geneList$locusName) {
  TimePlot(name)
  # print(TimePlot(name)$data$count)
}
dev.off()

#-------------------------combine DAM gene plots-------------------------
TimePlot <- function(x) {
  data <- plotCounts(dds, gene=x,intgroup=c("Stage","genotype"), returnData=TRUE)
  plot <-
    ggplot(data, aes(x=Stage, y=count, color=genotype, group=genotype)) +
    geom_point(size=2) + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
    scale_color_manual(values=c("#FF3333","#0033FF","#FF9999","#99CCFF"))
  #print(plot)
  #  return()
}


DAM1 = TimePlot("ppa018667m.g.v1.0")
DAM2 = TimePlot("ppb017585m.g.v1.0")
DAM3 = TimePlot("ppa010758m.g.v1.0")
DAM4 = TimePlot("ppa011123m.g.v1.0")
DAM5 = TimePlot("ppa010822m.g.v1.0")
DAM6 = TimePlot("ppa010714m.g.v1.0")
plot_DAM <- plot_grid(DAM1,DAM2, DAM3, DAM3,DAM5, DAM6, ncol = 3, nrow = 2)
plot_DAM + theme(text = element_text(size=8))
ggsave("DAMplots.png")
