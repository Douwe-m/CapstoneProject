library(Rsubread)

annotation_file <- "/local-fs/bachelor-students/2020-2021/Thema12/plant_gene_expression/iwgsc_refseqv2.1_annotation_200916_HC.gff3"
sam_file <- "/local-fs/bachelor-students/2020-2021/Thema12/plant_gene_expression/Aligned.out.sam"

output <- featureCounts(files = sam_file, 
                        annot.ext = annotation_file,
                        isGTFAnnotationFile = T,
                        GTF.featureType="mRNA",
                        GTF.attrType = "ID",
                        nthreads=64,
                        isPairedEnd=TRUE)
