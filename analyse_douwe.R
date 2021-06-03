library(Rsubread)

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
#Load the root count data
root_counts <- read.table("data/root_counts.csv", header = T, sep = ",", row.names = 1)


