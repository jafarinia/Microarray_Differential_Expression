# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
#+ Hossein Jafarinia modifications on R4.1.1 and Rstudio

#necessary libraries
library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(gplots)
library(reshape2)
library(plyr)
library(Biobase)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)

#Changing default VROOM_CONNECTION_SIZE so we can capture data
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 3)

#Making filing easy
curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub(paste0("/", sub("(.+)/", "", curD)), "", curD))

#load series and platform data from GEO
series = "GSE48558"
platform = "GPL6244"
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
# gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXX10XX10XX1X1",
               # "1X11X11XX1XX1XX1XX0X01XX0X0000X01X001X0010X010X010",
               # "XX10XX10XX1XXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               # "00000000000000000000")
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
smlWithName <- recode(sml, "1" = "AML", "0" = "Normal")
gset <- gset[ ,sel]


ex <- exprs(gset)

max(ex) 
min(ex)
#As we can see data or in log2

#Quality Control:
#Drawing Boxplot for quality control
gs <- factor(smlWithName)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE48558", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), legend = c("AML", "Normal"), bty="n")


corex = cor(ex)
pheatmap(corex, labels_row = smlWithName, labels_col = smlWithName, border_color = NA, main = "Smaples Correlation Heat Map") #As we can see Normal samples of different groups have highest correlation with each other and the lowest with AML samples and although AMS samples have some level of correlation with each other its not as much as normal samples due to the fact that cancer cells have some levels of difference between each other

#Dimensionality reduction
#PCA for genes
ex.scale <- t(scale(t(ex), scale = F)) #ex - mean(ex)
pc <- prcomp(ex.scale) #Finding PCAs
plot(pc, main = "PCA", xlab = "PCs")
#only PC1 and PC2 are sufficient and thats generally what our 2 dimensional page can show best
plot(pc$x[,1:2], main = "Genes PCA") #As we can see genes have a relatively meaningful distribution

#PCA for samples
pcr <- data.frame(pc$r[, 1:3], Group = smlWithName) #r means rotation we at least mean 3 pcs
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 3) + theme_bw() #As we can wee AML and normal cells are generally seperated well and we know from the examination summary says there are 5 types of Normal cells and we can see 5 clusters of normal cells in our plot which makes us conclude that the quality is good



#Differential Expression Analysis:
#assign samples to groups and set up design matrix
gset$group <- smlWithName
design <- model.matrix(~group + 0, gset) #finds pset
colnames(design) <- levels(gs) #design

#differential analysis main part by limma
fit <- lmFit(gset, design)  # fit linear model its a linear  model that fits a line to each model and based on the difference between line it say how much different they are

#set up contrasts of interest and recalculate model coefficients
cont.matrix <- makeContrasts(contrasts="AML-Normal", levels=design)
#cont.matrix and we see it wants to compare tumor and normal
fit2 <- contrasts.fit(fit, cont.matrix) #put the output in fit2

#compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01) #bayesian prior of 1% cancer related genes from biological studios
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf) # the function that calculated p.values and logFC "fdr" is our method for adjutment "B" is limmas test statistics and Inf for all genes

tT = subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC")) # a table of the columns we actually care about
write.csv(tT,"Results/tT.csv", row.names = FALSE)

#the genes that are expressed more in cancer samples(at least 2 times) their effect is meaningful
tT.Up.Gene <- subset(tT, adj.P.Val < 0.05 & logFC > 1)
upGenes = unique(tT.Up.Gene$Gene.symbol)
AML.Up.Genes.AllNames <- unique(as.character(strsplit2(upGenes, "///")))
write.table(AML.Up.Genes.AllNames, "Results/AMLUpGenes.txt", row.names = F, col.names = F, quote = F)

#the genes that are expressed less in cancer samples(at least 2 times) their effect is meaningful
tT.Down.Gene <- subset(tT, adj.P.Val < 0.05 & logFC < -1)
downGenes = unique(tT.Down.Gene$Gene.symbol)
AML.Down.Denes.AllNames <- unique(as.character(strsplit2(downGenes, "///")))
write.table(AML.Down.Denes.AllNames, "Results/AMLDownGenes.txt", row.names = F, col.names = F, quote = F)




#Extra work
# Multidimentional Scaling
library(magrittr)
library(dplyr)
library(ggpubr)
# Cmpute MDS
mds <- ex.scale %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
plot(mds) #as we can see its better