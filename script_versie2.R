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
rawdata <- read.table(file = "data/raw_counts.csv", sep = ",", header = TRUE)
rawdata[2:13] <- lapply(rawdata[2:13], as.numeric)

## Shows the first 5 lines of the table with rawdata. For the dimensions and structure of the data is used dim and str.
pander(rawdata)
# dim(rawdata)
# str(rawdata)

#Create a subset for every group
Leaf_Fe <- rawdata_visualizing[1:3]
Leaf_control <- rawdata_visualizing[4:6]
Root_Fe <- rawdata_visualizing[7:9]
Root_control <- rawdata_visualizing[10:12]

#3
## Visualizing the rawdata to a boxplot and a densityplot
## Boxplot

rawdata_log2 <- cbind(log2(rawdata[,2:13] + 1), rawdata[,1])
transformed_raw_data <- pivot_longer(rawdata_log2, cols = 1:12, names_to = "sample", values_to = "log2_value")

ggplot(transformed_raw_data, aes(x = sample, y = log2_value))+ 
  geom_boxplot(fill = "steelblue") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle  = 45, hjust = 1)) + 
  labs(x = "Sample", y = "log2(counts)")

## Density plot
myColors <- hue_pal()(4)
plotDensity(log2(rawdata[,-1] + 0.1), col=rep(myColors, each=3),
            lty=c(1:ncol(rawdata[,-1])), xlab='Log2(count) of the raw data',
            ylab = 'Density', main='Density plot')
legend('topright', names(rawdata[,-1]), lty=c(1:ncol(rawdata[,-1])),
       col=rep(myColors, each=3))

## Barplot
col_sums <- as.data.frame(apply(rawdata[2:13], 2, sum))
col_sums$group <- row.names(col_sums)
rownames(col_sums) <- 1:12
colnames(col_sums) <- c("sum", "group")

ggplot(col_sums, aes(x = group, y = sum/1e6)) + 
  geom_bar(stat = "identity", fill = "steelblue") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle  = 45, hjust = 1)) + 
  labs(x = "Sample", y = "Total read count in millions")

## Heatmap
ddsMat <- DESeqDataSetFromMatrix(countData = rawdata_visualizing,
                                 colData = data.frame(samples = names(rawdata_visualizing)),
                                 design = ~ 1)
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

groups <- factor(rep(1:4, each=3), 
                labels = c("leaf_treatment", "leaf_control", "root_treatment", "root_control"))


names(mdsPoisData) <- c('x_coord', 'y_coord')

Legend <- factor(rep(1:4, each=3), 
                 labels = c("Leaf_Fe", "Leaf_control", "Root_Fe", "Root_control"))
coldata <- names(rawdata_visualizing)

ggplot(mdsPoisData, aes(x_coord, y_coord, color = groups, label = coldata)) + 
  geom_text(size = 4) +
  ggtitle('Multi Dimensional Scaling') +
  labs(x = "Poisson Distance", y = "Poisson Distance") +
  theme_light()

#4
## Pre Processing data for DEG analysis
counts.fpm <- log2( fpm(ddsMat, robust = TRUE) + 1 )

#DESeq2

dds <- DESeqDataSetFromMatrix(countData = rawdata,
                              colData = data.frame(condition=condition),
                              design = ~ 0+condition)
results_dds <- DESeq(dds)

resultsNames(results_dds)
res <- results(results_dds, contrast = c("condition", "Leaf_Fe", 
                                         "Leaf_Control"))
res1 <- results(results_dds, contrast = c("condition", "Root_Fe", 
                                          "Root_Control"))
summary(res)

## Volcano

# Load packages
library(gridExtra)
# Determine the p-value threshold and log fold change 
deseq.results <- cbind(rld, as.data.frame(res@listData))
pval_threshold <- 0.05
logfc_threshold <- 1
deseq.threshold <- as.factor(abs(deseq.results$log2FoldChange) >= logfc_threshold & 
                               deseq.results$padj < pval_threshold)
# Plot the data to a volcano plot
g = ggplot(data=deseq.results, 
           aes(x=log2FoldChange, y=-log10(padj), 
               colour=deseq.threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  theme_bw() + theme(legend.position="none") +
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle("Leaf_Fe vs Leaf_Control DEGs according to DESeq2 with adjusted 
          p-value <= 0.05 and absolute FC >= 2") +
  theme(plot.title = element_text(size = 8)) +
  geom_text_repel(aes(label=ifelse(padj < 0.1e-6 & abs(log2FoldChange) 
                                   >= logfc_threshold, row.names(deseq.results), 
                                   '')))
# Make the data compatible for sample 2
deseq.results1 <- cbind(rld, as.data.frame(res1@listData))
deseq.threshold1 <- as.factor(abs(deseq.results1$log2FoldChange) 
                              >= logfc_threshold & deseq.results1$padj < pval_threshold)

g1 = ggplot(data=deseq.results1, 
            aes(x=log2FoldChange, y=-log10(padj),
                colour=deseq.threshold1)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  theme_bw() + theme(legend.position="none") +
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle("Root_Fe vs Root_Control DEGs according to DESeq2 with FDR <= 0.05 and absolute FC >= 2") +
  theme(plot.title = element_text(size = 8)) +
  geom_text_repel(aes(label=ifelse(padj < 0.1e-6 & abs(log2FoldChange) 
                                   >= logfc_threshold,
                                   row.names(deseq.results1), '')))

grid.arrange(g, g1, nrow = 1) 

# Venn Diagram

DESeq.degs <- row.names(deseq.results[which(deseq.threshold == TRUE),])
DESeq.degs1 <- row.names(deseq.results1[which(deseq.threshold1 == TRUE),])

library(VennDiagram)

deg.intersect = length(intersect(DESeq.degs, DESeq.degs1))
deg.venn <- list('intersect' = deg.intersect,
                 'Leaf' = length(DESeq.degs),
                 'Root' = length(DESeq.degs1))


venn.plot <- draw.pairwise.venn(deg.venn$Leaf, deg.venn$Root, deg.venn$intersect,
                                category = c("Leaf", "Root"), scaled = F,
                                fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                cat.pos = c(0, 0))
grid.draw(venn.plot)

# Actually plot the plot
grid.draw(venn.plot)

library(dplyr)

deseq.results %>% arrange(padj) %>% arrange(desc(log2FoldChange)) %>% head(10)