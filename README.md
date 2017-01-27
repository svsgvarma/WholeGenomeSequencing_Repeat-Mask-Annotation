# WholeGenomeSequencing_Repeat-Mask-Annotation

-Whole Genome Sequcing (WGS) of individual sample VCF file or multi-sample VCF file can be recruted for annotation pipleline, this allows to check for varaints falling with in the low-complexity regions

-Script to filter out false negatives and true negatives of variants with in the Repeat-Masking(RM) regions, RM data was obtained from UCSC browser

-Checks and retain reliable variants which are falling within concordant regions of public resources dowloaded for ~200 WGS data

-Check against local 10 WGS datasets for futher filtering step 

-Finally another script to add gene expression FPKM scores from Gene-Epression-BodyMap2.0 public resource for 12 different tissues, annotations for each gene can be added
