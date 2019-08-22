setwd("~/Desktop/Jiali/UTK/Peach/")
#install.packages("Rmisc")
library("Rmisc")
if(!require(lsmeans)){install.packages("lsmeans")}
install.packages("multcompView")
library(lsmeans)
library(multcompView)
library("ggplot2")

budData <- read.csv("bud size.csv", header = T, row.names = 1)
budData$Timepoint <- as.factor(budData$Timepoint)
width <- summarySE(budData, measurevar = "Width", na.rm = TRUE,
                   groupvars = c("Genotype","Timepoint"))
width$Genotype <- factor(width$Genotype, levels = c("A209","A340","A318","A323"))

length <- summarySE(budData, measurevar = "Length", na.rm = TRUE,
                   groupvars = c("Genotype","Timepoint"))
length$Timepoint <- as.factor(length$Timepoint)
length$Genotype <- factor(length$Genotype, levels = c("A209","A340","A318","A323"))

carpel <- summarySE(budData, measurevar = "Carpel", na.rm = TRUE,
                   groupvars = c("Genotype","Timepoint"))
carpel$Timepoint <- as.factor(carpel$Timepoint)
carpel$Genotype <- factor(carpel$Genotype, levels = c("A209","A340","A318","A323"))


aov(Width ~ Genotype+Timepoint+Genotype:Timepoint, data = budData)
p<- ggplot(width, aes(x=Timepoint, y=Width, fill=Genotype)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Width-sd, ymax=Width+sd), width=.2,
                position=position_dodge(.9)) +scale_fill_grey(breaks=c("A209","A340","A318","A323"))
p+labs(title="Peach Bud width", x="Chilling hours (hrs)", y = "Width (mm)") + theme_classic()


model = lm(Carpel ~ Genotype + Timepoint + Genotype:Timepoint,
           data = budData)
anova(model)
PairCom = lsmeans(model,
                   ~ Genotype:Timepoint)
CLD = cld(PairCom,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")

CLD
CLD$Genotype <- factor(CLD$Genotype, levels = c("A209","A340","A318","A323"))
CLD$.group=gsub(" ", "", CLD$.group)

p<- ggplot(CLD, aes(x=Timepoint, y=lsmean, fill=Genotype)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=.2,
                position=position_dodge(.9)) +scale_fill_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF"))+
  geom_text(aes(label = .group, y = lsmean + 0.2),
            position = position_dodge(0.9),
            vjust = 0,
            color   = "black")
p+labs(title="Peach Bud Carpel Length", x="Chilling hours (hrs)", y = "Length (mm)") + theme_classic(base_size = 16)

