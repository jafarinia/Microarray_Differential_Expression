# library(limma)
# x <- matrix(runif(30, 0, 15), nrow = 10)
# rownames(x) <- paste0("G", 1:10)
# colnames(x) <- paste0("Sample", 1:3)
# y <- normalizeQuantiles(x) # quantile normalization
# pc <- prcomp(y) # find all principal component analyzises
# plot(pc)
# plot(pc$x[,1:2]) #plot it on pcs



library(limma)
library(pheatmap)
library(ggplot2)
x <- matrix(runif(80, 0, 15), nrow = 10)
rownames(x) <- paste0("G", 1:nrow(x))
colnames(x) <- c(paste0("Normal", 1:4), paste0("Tumor", 1:4))
x[1:3, 1:4] <- x[1:3, 1:4] + rnorm(12, mean = 12) # artificially making first 3 genes in normals have higher expression 
x[4:6, 5:8] <- x[4:6, 5:8] + rnorm(12, mean = 12) # artificially making second 3 genes in in cancers have higher expression 
pheatmap(x)
y <- normalizeQuantiles(x) # quantile normalization
pheatmap(y)
pc <- prcomp(y) # find all principal component analyzises
plot(pc)
plot(pc$x[,1:2]) #plot it on pcs 
pcr <- data.frame(pc$r)#We did this because the number of samples is lower than then number of genes so pcr shows samples
pcr$Group <- c(rep("Normal", 4), rep("Tumor", 4))
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 4) + theme_bw()


dif <- data.frame(Normal = rowMeans(y[,1:4]), Tumor = rowMeans(y[,5:8])) #Averaging over first 4 and second 4
#For example rowMeans(x) gives is average over each column because defualt is average is pver each column
dif$logFC <- dif$Normal - dif$Tumor
t.test(y[1,1:4], y[1, 5:8]) # student t test on gene 1
grades <- c(15, 25, 20, 30, 20, 5, 35, 16, 11, 35, 20.5, 12, 30, 30, 27, 30, 20, 9, 20, 9, 5, 25, 25, 10, 20, 25, 15, 7,           20, 20, 30, 5, 12, 15, 18, 14, 9, 8)
maxgrades <- c(30, 35, 35, 30, 20, 20.5, 30, 20, 18)
t.test(maxgrades[1:7], maxgrades[8:9])