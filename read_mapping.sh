#Stap 1: downloaden SRA bestand voor 1 sample
#Stap 2: trimmen en filteren
#Stap 3: mappen, output .BAM
#Stap 4: .fastq verwijderen	

sra_accession = SRR13114670

#Download de .fastq bestanden van het sample
fasterq-dump $sra-accession

#Trimmen en filteren
java -jar /local-fs/bachelor-students/2020-2021/Thema12/plant_gene_expression/Trimmomatic-0.39/trimmomatic-0.39.jar 
	PE $sra_accession_1.fastq sra_accession_2.fastq 
	trimmed_$sra_accession_1_p.fastq trimmed_$sra_accession_1_u.fastq 
	trimmed_$sra_accession_2_p.fastq trimmed_$sra_accession_2_u.fastq 
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:1 MINLEN:36 SLIDINGWINDOW:4:20 AVGQUAL:20
