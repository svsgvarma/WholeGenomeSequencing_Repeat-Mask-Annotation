#!/bin/bash



###-----
#./runANNOTAIONS_all-samples.bash > ./work_WGS/runANNOTAIONS_all-samples.bash.log 2>&1

#####------Calling correct protein CODON -------
#Annotation_Coding-region-check.py INPUT1=FULLPATH+Inputfile(VCF) INPUT1=Outputfile1 INPUT3=Outputfile2 INPUT4=Working-Directory-fullpath/

##----DCM: WGS and WES----

WORKDIR="/media/varma/Maxtor/ANNOTATION_Scripts/"

WORK_ANN1="/media/varma/SAMSUNG/HCM_WGS_ANNOTATION_ANNOVAR/"
WORKDIR="/media/varma/Maxtor/ANNOTATION_Scripts/"

python "$WORKDIR"Annotation_Coding-region-check.py /media/varma/My_Passport_2/Cardiomyopathy_WES_BACKUP/CMP_WES_ANNOTATION_ANNOVAR/samp12345_withoutgel_1720017671182881920819237_variants_GATK_step2_annovar.hg19_multianno.vcf samp12345_withoutgel_1720017671182881920819237_Protein-CODON_changes.txt $WORK_ANN1

python "$WORKDIR"Annotation_Coding-region-check.py /media/varma/My_Passport_2/Cardiomyopathy_WES_BACKUP/CMP_WES_ANNOTATION_ANNOVAR/samp12345_withgel_1720017671182881920819237_variants_GATK_step2_annovar.hg19_multianno.vcf samp12345_withgel_1720017671182881920819237_Protein-CODON_changes.txt $WORK_ANN1

python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN1"HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf HCM_B00H7EW-B00H7EX-B00H7EY_Protein-CODON_changes.txt $WORK_ANN1

##----ADCA: WGS ----

#WORK_ANN2="/media/varma/Maxtor/ADCA_WGS_ANNOTATION_ANNOVAR/"
#WORKDIR="/media/varma/Maxtor/ANNOTATION_Scripts/"
#python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN2"ADCA_009-028-080_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf ADCA_009-028-080_Protein-CODON_changes.txt $WORK_ANN2


##----MALRARWGSTEV:WGS----

WORK_ANN3="/media/varma/Maxtor/MALRARWGSTEV_ANNOTATION/"

python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN3"AAD-704_B00HP0S-B00HP0T_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf AAD-704_B00HP0S-B00HP0T_Protein-CODON_changes.txt $WORK_ANN3

python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN3"SAL-REZ-394_B00HP0Q-B00HP0R_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf SAL-REZ-394_B00HP0Q-B00HP0R_Protein-CODON_changes.txt $WORK_ANN3

python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN3"SAL-VAN-379_B00HP0O-B00HP0P_chr1-MT_variants_GATK_step2_annovar.hg19_multianno.vcf SAL-VAN-379_B00HP0O-B00HP0P_Protein-CODON_changes.txt $WORK_ANN3

echo "Done running script... Protein CODAN changes...."

###-----
##----DVT-6samp:WES----

WORK_ANN3="/media/varma/My_Passport_2/DVT_WES_ANNOTATION_ANNOVAR/"
WORKDIR="/media/varma/Maxtor/ANNOTATION_Scripts/"
#python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN3"samp123456_FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6_variants_GATK_step2_annovar.hg19_multianno.vcf samp123456_FSCD1FSCD2FSCD3FSCD4FSCD5FSCD6_Protein-CODON_changes.txt $WORK_ANN3

#python "$WORKDIR"Annotation_Coding-region-check.py "$WORK_ANN3"samp789101112_FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12_variants_GATK_step2_annovar.hg19_multianno.vcf samp789101112_FSCD7FSCD8FSCD9FSCD10FSCD11FSCD12_Protein-CODON_changes.txt $WORK_ANN3


########################


