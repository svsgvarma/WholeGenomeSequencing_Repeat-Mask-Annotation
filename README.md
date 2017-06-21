# WholeGenomeSequencing_Repeat-Mask-Annotation

-Whole Genome Sequcing (WGS) of individual sample VCF file or multi-sample VCF file can be recruted for annotation pipleline, this allows to check for varaints falling with in the low-complexity regions

-Script to filter out false negatives and true negatives of variants with in the Repeat-Masking(RM) regions, RM data was obtained from UCSC browser

-Checks and retain reliable variants which are falling within concordant regions of public resources dowloaded for ~200 WGS data

-Check against local 10 WGS datasets for futher filtering step 

-Finally another script to add gene expression FPKM scores from Gene-Epression-BodyMap2.0 public resource for 12 different tissues, annotations for each gene can be added

# Scripts for annotation of WGS

# 1. Annotation pipeline:

After variant calling from GATK, I have annotated with annovar and resulting vcf file was generated (File extension: *_variants_GATK_step2_annovar.hg19_multianno.vcf). To improve further, I have performed series of additional annotation steps as follows. Diagram illustrates the full annotation pipeline in the Figure Supplementary S2.

1) Parsing raw data file:
Initially, I have parsed raw annovar output file (File extension: *_variants_GATK_step2_annovar.hg19_multianno.txt) and retrieved only necessary columns, and generated a new file (File extension: *_Variants.txt). This contains the whole variants (file size ~2.2 GB, 7058947 variants).

Since, DCM disease is dominant type and there are no control samples, I have selected only heterozygous (0/1; 1/0), and  generated just heterozygous filtered file (File extension: *_Variants_Het.txt).

2) Applied Basic filters: 
I made three initial standard filters to reduce the number as below.
2.1) Along with heterozygous filter(Het=Affected Hom=Unaffected)
2.2) QUAL (VQSR) == PASS filter (https://software.broadinstitute.org/gatk/guide/article?id=39). For this two versions of files were made, one with PASS filter another one without PASS filter. 
2.3) MAF<5% (~13MB = 39543 variants)
For the next step, I have marked and added below columns without removing any variants from the above file (~13MB = 39543 variants). 
( File extension: *_Variants_Het_MAF5%.txt)
( File extension: *_Variants_Het_PASS_MAF5%.txt)

3) Prediction Scores:
There are several prediction scores have been generated using annovar annotation, however, for WGS newly published scoring systems, especially complete CADD datasets for indels are not available in annovar. Therefore, I have downloaded all datasets locally and annotated with CADD, ReMM, LISIGHT and DANN scores (File extension:*_Variants_Het_PASS_MAF5%_CRLD_RE-RG-LDB.txt). I have also made another file with these scores without applying above basic filters 
(File extension:*_Variants_Het_CRLD.txt).

About prediction score thresholds:
For ReMMM scores, the authors recommend "probably damaging" if the score is above 0.9 and close to 1.0. Except ReMM scores, I do not see any recommended thresholds for other prediction scores.

They are usually in the same range between +0 to +1, score close to 0 is not good (common; non-damaging), score close to 1 is a good (rare; probably damaging), except few like Eigen, FATHMM  and CADD scores ranges from -10 to +10. Maybe, it is not good idea to filter on Eigen or below scores because each scoring system has pros and cons, even sometimes good variant might have -ve scores in one and +ve in another. please check example below. 

3.1) Eigen: Range= -10; 0 to +10 (I think, eigen value smaller than -2.0 or even smaller might not be a good one)
3.2) CADD1.3_PHRED: Range= +0 to +35 (even grater than this, and higher score is better)
3.3) CADD1.3_RawScore:  Range= -10 to +10
3.4) FATHMM_noncoding: Range= -99; +0 to +1
3.5) FATHMM_coding: Range= -99; +0 to +1
3.6) ReMM_Scores: Range= 0-1 (over 0.9 is probably damaging)
3.7) LINSIGHT_Scores: Range= 0-1
3.8) DANN_Scores: Range= 0-1

Example: Eigen:-0.5518 CADD1.3_RawScore:0.640551    CADD1.3_PHRED:8.433    ReMM_Scores:0.881     LINSIGHT_Scores:.    DANN_Scores:0.9855

4) Reliable Genome (from 219 WGS control samples):
Several persisting False positives and negatives should be removed, one possible solution would be to remove variants which are falling with in Repeat Masking (RM) regions and discordant regions, this can be achievable by some recently published data which is ReliableGenome [15], generated based on 219 WGS datasets. I have passed variants to this data and marked them as below.
4.1) RG = ReliableGenome, marked as RG. (Niko Popitsch, et,al)
4.2) RG = ReliableGenome_min-LCR-min-HD, marked as RG_min-LCR-HD. (Niko Popitsch, et,al)
4.3) RG = ReliableGenome_min-UM75, marked as RG_min-UM75. (Niko Popitsch, et,al)
( File extension: *_Variants_Het_PASS_MAF5%_CRLD_RE-RG-LDB.txt)

5) Repeat Elements:
In addition, I assume there must be more Repeat Elements (RE) data available in UCSC browser, thereby I have downloaded another set of RE data and developed in house script to mark these regions too.
RE = Repeat Elements downloaded file from UCSC browser. Marked as NotIn_RE (variant is Not present In the Repeat Elements) or In_RE (variant is Present In the Repeat Elements). 
( File extension: *_Variants_Het_PASS_MAF5%_CRLD_RE-RG-LDB.txt)

6) Local WGS data search:
As a control dataset, I have searched each variant position against already called variants from other disease datasets. Same pipeline, with HG19 was used to call variants in all other control samples. 
LDB = Local Data Base search was performed on other 63 WES (SCA-23F, DVT-3F, HCM-1F, DCM-1F, DD-LD-16F) and 22 WGS samples (of 3-DCM, 5-HHT1,5-HHT2). 
Example: marked as NotIn_LDB-WGS and In_LDB-WGS-H2 etc. 
(File extension: *_Variants_Het_PASS_MAF5%_CRLD_RE-RG-LDB.txt)

7) Gene Expression data:
Gene Expression FPKM scores was applied for all variants from two databases as below.
7.1) Gene-Epression-BodyMap2.0 = #Gene-Epression-BodyMap2.0_16-tissues, following the gene name column 16 columns of FPKMs scores (http://www.cureffi.org/2013/07/11/tissue-specific-gene-expression-data-based-on-human-bodymap-2-0/) were added for each tissue type.

7.2) Gene-Epression-GTEx = #Gene-Epression-GTEx_53-tissues, included expression RPKM scores for 53 tissues (http://gtexportal.org/home/tissueSummaryPage).
(File extension: *_CRLD_RE-RG-LDB_GExp_REDIportal_GnC_UCSC.txt)

8) REDIportal:
REDIpotal [16] data  was downloaded and applied with 24 columns after Gene expression columns, there is also a user friendly web-browser (http://srv00.recas.ba.infn.it/atlas/index.html) available for this data.
(File extension: *_CRLD_RE-RG-LDB_GExp_REDIportal_GnC_UCSC.txt)

9) Variant count on each gene:
Here, I made another column Gene Count (#Gene based Variants count), showing the count of variants falling on each unique gene. Example shown below, CHD5 gene has 4 variants at this filtering stage (MAF 5%). Basic idea is that, one can tract the number of variants falling on a gene level, even after further filtering steps. Also, this could be useful to check for two rare heterozygous variants with in the same gene. 
(File extension: *_CRLD_RE-RG-LDB_GExp_REDIportal_GnC_UCSC.txt)
#Gene based Variants count
#Gene_Names
4
"CHD5"
4
"CHD5"
4
"CHD5"
4
"CHD5"
2
"RPL22"
2
"RPL22"

10) UCSC format file:
UCSC supported input BED format, added 9 columns at the beginning of the file. After user filtering, just only 9 columns can be copied and added to the UCSC browser as an input bed format. Proper header should be added to these columns as example shown below to support for the UCSC browser required formats (https://genome.ucsc.edu/goldenpath/help/customTrack.html).

Example_header: "track name='DCM_WGS' description='DCM_three-samples' visibility=2 itemRgb=On". Diagram shown in the Figure Supplementary S3.
(File extension: *_CRLD_RE-RG-LDB_GExp_REDIportal_GnC_UCSC.txt)

10.1. RGB (Read Green Blue) color indication code as follows:
####------Exonic---------
'255,0,0' == ('exonic' , ExAC_ALL == "." , 1000GP == "."  and GnomAD == "." )
'255,0,55' == ('exonic')
####------Non-Exonic---------
'0,0,255' == ('non-exonic' ,  1000GP == "." and GnomAD == ".")
'0,55,255' == ('non-exonic')

11) Null or missing genotype: 
Included null calls or missing genotype variants (./.). GATK caller fail to call these variants because of their bad coverage or lack of reads. I made a separate file including these variants.
(File extension: *_Variants_Het+N_MAF5%.txt)

12) Protein codon changes: 
Annotation of paired or even continuous SNVs (including synonymous and non-synonymous variants) located within the same genetic code requires attention. This idea conflicts with the annotations observed by traditional annotation software widely used now a days (ex: ANNOVAR tool). While looking at the combined effect within the framework of genetic code, we have new coding amino acid predicted by using in-house python scripts.
(File extension: *_Protein-CODON_changes.txt)
