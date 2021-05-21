# CapstoneProject

## Introductie

Een tekort aan ijzer kan de gewasopbrengst sterk beïnvloeden. Daarnaast speelt ijzer een belangrijke rol in het metabolisme en is het nodig voor plantengroei. IJzer fungeert onder andere als cofactor van veel enzymen. Verder is het betrokken bij de elektronentransportketen en fotosynthese. In eerdere onderzoeken is al onderzoek gedaan naar ijzertekorten in verschillende planten, echter is het achterliggend moleculaire mechanisme betreffende ijzertekorten niet bekend voor tarwe. Met behulp van RNA-seq zal getracht worden dit mechanisme te ontrafelen. [(Wang et al., 2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7550799/)

De gebruikte cultivar is Triticum aestivum cv. Bobwhite S26 s. Zaden werden in gelijke omstandigheden gekiemd met voldoende ijzer aanwezig. Na 7 dagen werd de helft van de zaden overgebracht naar een ijzerarme omgeving. Vervolgens werden de planten gegroeid gedurende 90 dagen. Vervolgens werden van de wortels en bladeren van de planten in de twee condities RNA geïsoleerd. Per sample werden drie biologische replicaten genomen. Na een zuiveringsstap werden de sequencing libraries gemaakt. Na de library preparation werden de samples met behulp van de Illumina HiSeq gesequenced. In totaal werden 12 samples gesequenced. [(Wang et al., 2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7550799/)

 ## Stappen uit het artikel
 
- QC: fastqc
- Trimmen: Trimmomatic-0.36 (cutting Illumina TruSeq adapters from the reads)
- Mappen: STAR 2.7.3a (allowing maximum 5 bp mismatches, WGSC wheat genome assembly RefSeq v1.0)
- featureCounts v1.6.4
- DEseq2 1.26.0

## To Do
 - [X] SRA data downloaden
 - [X] QC
 - [X] [Genoom indexeren](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
 - [X] [Mappen met STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
 - [ ] Homeoroq?
 - [ ] featureCounts
