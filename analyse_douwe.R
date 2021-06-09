#Load the required packages
library(Rsubread)
library(DESeq2)
library(EnhancedVolcano)

#----------------------------------------------------------------------
#Path to the annotation file
annotation_file <- "/local-fs/bachelor-students/2020-2021/Thema12/plant_gene_expression/iwgsc_refseqv2.1_annotation_200916_HC.gff3"

#Paths to the .BAM files
bam_files <- paste0("/local-fs/bachelor-students/2020-2021/Thema12/plant_gene_expression/raw_data/sample_SRR131146", 
                    70:75, 
                    "/Aligned.sortedByCoord.out.bam")

#Run featureCounts
featureCounts_output <- featureCounts(files = bam_files, 
                                      annot.ext = annotation_file,
                                      isGTFAnnotationFile = T,
                                      GTF.featureType="mRNA",
                                      GTF.attrType = "ID",
                                      nthreads=64,
                                      isPairedEnd=TRUE)

#Save the count data in a separate dataframe and add column names
count_data <- as.data.frame(featureCounts_output$counts)
colnames(count_data) <- c("root_treatment_1", "root_treatment_2", "root_treatment_3", 
                          "root_control_1", "root_control_2", "root_control_3")

#Save the count data as a csv file
write.csv(count_data, "data/root_counts.csv")

#----------------------------------------------------------------------
#Load root count data
root_counts <- read.table("data/root_counts.csv", header = T, sep = ",", row.names = 1)

#sample information table
coldata <- data.frame(Treatment = as.factor(rep(c("Low_Fe", "Control"), each = 3)), 
                      row.names = colnames(root_counts))

#construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = root_counts,
                              colData = coldata,
                              design = ~ Treatment)

#Differential expression analysis
results_dds <- DESeq(dds)
res  <- results(results_dds)

#Calculate the amount of genes that are up/down regulated
deseq.threshold <- as.factor(abs(res$log2FoldChange) >= 1 & res$padj < 0.05)
deseq.degs <- row.names(res[which(deseq.threshold == TRUE),])
df <- data.frame(up = sum(res[deseq.degs,]$log2FoldChange > 1, na.rm = T),
                 down = sum(res[deseq.degs,]$log2FoldChange < -1, na.rm = T),
                 total = length(deseq.degs),
                 row.names = "Root treatment vs root control")
knitr::kable(df)


#Visualise the results
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue', 
                pCutoff = 1e-05,
                FCcutoff = 1,
                labSize = 0,
                legendPosition = "none",
                subtitle = "root treatment vs root control", 
                gridlines.major = F, 
                gridlines.minor = F)