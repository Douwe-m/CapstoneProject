# load libraries
library(ggplot2)
library(pander)
library(pheatmap)
library(affy)
library(scales)
library(DESeq2)
library(PoiClaClu)
library(tidyr)

#2

## loading CSV with raw counts in R, and make the integer columns numeric.
rawdata <- read.table(file = "raw_counts.csv", sep = ",", header = TRUE)
rawdata[2:13] <- lapply(rawdata[2:13], as.numeric)

## Shows the first 5 lines of the table with rawdata. For the dimensions and structure of the data is used dim and str.
pander(rawdata)
dim(rawdata)
str(rawdata)

## Split the data for the different samples (leaf, leaf_control, root, root_control)
print(names(rawdata))
row.names(rawdata) <- rawdata$Gene

rawdata_visualizing <- rawdata[-1]
names(rawdata_visualizing)

Leaf_Fe <- rawdata_visualizing[1:3]
Leaf_control <- rawdata_visualizing[4:6]
Root_Fe <- rawdata_visualizing[7:9]
Root_control <- rawdata_visualizing[10:12]

#3
## Visualizing the rawdata to a boxplot and a densityplot
summary(rawdata_visualizing)


## Boxplot

rawdata_log2 <- log2(rawdata[,2:13] + 1)
rawdata_gene_names <- rownames(rawdata)
rawdata_boxplot <- cbind(rawdata_gene_names, rawdata_log2)

rawdata_boxplot %>% 
  gather(Leaf_Fe_1, rawdata_gene_names) %>% 
  ggplot(aes(Leaf_Fe_1, rawdata_gene_names)) + 
  geom_boxplot()

## Density plot
myColors <- hue_pal()(4)
plotDensity(log2(rawdata_visualizing + 0.1), col=rep(myColors, each=3),
            lty=c(1:ncol(rawdata_visualizing)), xlab='Log2 count of the raw data',
            ylab = 'Density', main='Density plot')
legend('topright', names(rawdata_visualizing), lty=c(1:ncol(rawdata_visualizing)),
       col=rep(myColors, each=3))

barplot(colSums(rawdata_visualizing) / 1e6)

ddsMat <- DESeqDataSetFromMatrix(countData = rawdata_visualizing,
                                 colData = data.frame(samples = names(rawdata_visualizing)),
                                 design = ~ 1)

## Heatmap
rld.dds <- vst(ddsMat)
rld <- assay(rld.dds)
sampledists <- dist( t( rld ))
sampleDistMatrix <- as.matrix(sampledists)

annotation <- data.frame(Treatment = factor(rep(rep(1:2, each = 3), 2),
                                            labels = c("Fe", "Control")),
                         Tissue = factor(rep(1:2, each = 6), 
                         labels = c("Leaf", "Root")))

rownames(annotation) <- names(rawdata_visualizing)

pheatmap(sampleDistMatrix, show_colnames = FALSE,
         annotation_col = annotation,
         clustering_distance_rows = sampledists,
         clustering_distance_cols = sampledists,
         main = "Headmap")

## MDS plot
dds <- assay(ddsMat)
poisd <- PoissonDistance( t(dds) )
samplePoisDistMatrix <- as.matrix(poisd$dd)
mdsPoisData <- data.frame( cmdscale(samplePoisDistMatrix) )

names(mdsPoisData) <- c('x_coord', 'y_coord')

Legend <- factor(rep(1:4, each=3), 
                 labels = c("Leaf_Fe", "Leaf_control", "Root_Fe", "Root_control"))
coldata <- names(rawdata_visualizing)

ggplot(mdsPoisData, aes(x_coord, y_coord, color = groups, label = coldata)) + 
  geom_text(size = 4) +
  ggtitle('Multi Dimensional Scaling') +
  labs(x = "Poisson Distance", y = "Poisson Distance") +
  theme_bw() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

#4
## Pre Processing data for DEG analysis

counts.fpm <- log2( fpm(ddsMat, robust = TRUE) + 1 )
















