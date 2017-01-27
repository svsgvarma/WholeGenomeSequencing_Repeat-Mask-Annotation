#!/usr/bin/env python


import sys
import subprocess


#############------------

"""

WholeGenomeSequencing_Repeat-Mask-Annotation

-Whole Genome Sequcing (WGS) of individual sample VCF file or multi-sample VCF file can be recruted for annotation pipleline, this allows to check for varaints falling with in the low-complexity regions


INPUT_FILES:

python checkfor_UCSC_RG-LCR-HD-UM75_LocalDB.py /home/varma/proj/HCM_WGS_ANNOTATION_ANNOVAR_ReMM/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_GATK_ANNOVAR_het_PASS-5%.txt /media/varma/My_Passport_2/Scripts-for-LocalDBsearch/OUT_DATA/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_GATK_ANNOVAR_het_PASS-5%_notinUCSC-inRG.txt /media/varma/My_Passport_2/Scripts-for-LocalDBsearch/OUT_DATA/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_GATK_ANNOVAR_het_PASS-5%_notinUCSC-inRG_notinLocalDB.txt

"""

class fileHandler:
	def __init__(self):
		self.data = []
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):
	def __init__(self):
		self.data = []
		from collections import defaultdict
		self.ident_ranges_UCSC = defaultdict(list)
		self.ident_ranges_RG = defaultdict(list)
		self.ident_ranges_RG_LCR_HD = defaultdict(list)
		self.ident_ranges_RG_UM75 = defaultdict(list)
		
	def Lsearch_UCSC(self,readfl1,writefl3):
		"""
		Calling UCSC localsearch method:
		-Script to filter out false negatives and true negatives of variants with in the Repeat-Masking(RM) regions, RM data was obtained from UCSC browser
		"""

		self.ident_ranges = defaultdict(list)
		with open("/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/RG-ReliableGenome/Repeat-masking/rmsk-pos.txt",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split()
				start, end = int(start), int(end)
				self.ident_ranges[ident[3:]].append((start, end))
		with open(readfl1,'r') as f1, open(writefl3,'w') as output:
			first_line = f1.readline()
			output.write(first_line)
			for line in f1:		
				ident, value = line.strip().split()[0:2]
				value = int(value)
				if not any(start <= value <= end for start, end in self.ident_ranges[ident]):
					output.write(line)
		return None

	def Lsearch_UCSC_RG_LCR_HD_UM75(self,readfl1,writefl3):
		"""
		Calling RG localsearch method
		-Checks and retain reliable variants which are falling within concordant regions of public resources dowloaded for ~200 WGS data

		"""
		DIR="/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/RG-ReliableGenome/wtchg-rg-master/data/public/"
		with open("/media/varma/My_Passport_2/Scripts-for-LocalDBsearch/RG-ReliableGenome/Repeat-masking/rmsk-pos.txt",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split()
				start, end = int(start), int(end)
				self.ident_ranges_UCSC[ident[3:]].append((start, end))
		with open(DIR+"20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG[ident].append((start, end))
		with open(DIR+"20160422_RG-win1000-score1_-3-RELIABLE-above0.5.bed.gz-min-LCR-min-HD.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG_LCR_HD[ident].append((start, end))
		with open(DIR+"20160825_RG-win1000-score1_-3-RELIABLE-above0.5.-min-UM75.bed",'r') as f2:
			for row in f2:
				ident, start, end = row.strip().split("\t")[0:3]
				start, end = int(start), int(end)
				self.ident_ranges_RG_UM75[ident].append((start, end))

		with open(readfl1,'r') as f1, open(writefl3,'w') as output:
			first_line = f1.readline()
			output.write(first_line)
			for line in f1:
				ident, value = line.strip().split()[0:2]
				value = int(value)
				if not any(start <= value <= end for start, end in self.ident_ranges_UCSC[ident]):
					True
					if any(start <= value <= end for start, end in self.ident_ranges_RG[ident]):
						True
						if any(start <= value <= end for start, end in self.ident_ranges_RG_LCR_HD[ident]):
							True
							if any(start <= value <= end for start, end in self.ident_ranges_RG_UM75[ident]):
								output.write(line)

		print "Done seach for UCSC_RG_LCR_HD_UM75...."

		return None

	def Lsearch_WGS_WES(self,rfl1,wfl2):
			"""
			Calling localsearch of  WGS and WES method
			-Final check against local 10 WGS and ~ 90 WES individual datasets for futher filtering step

			"""
			infl1 = self.open_file(rfl1) 
			infl_nhr = infl1[1:]
			file2 = self.write_file(wfl2)
			file2.write(str(infl1[0]))
			totlns= len(infl_nhr)
			count1 = 0;
			while (count1 < totlns):
				infll = infl_nhr[count1].strip().split()
				#grep -w '1'$'\t''15118'
				try:
					True
					cmdFls1 = "for i in /media/varma/My_Passport_2/*_WES_VCF/*_SNPnINDELcalls.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
					grepout1 =  subprocess.check_output(cmdFls1, shell=True)
				except:
					True
					#print "not present in 39_DB"
					try:
						True
						cmdFls2 = "for i in /media/varma/My_Passport_1/Backup_Varma/RenduOsler1/work_WGS/39_12345_*_variants_Varscan_SNP.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
						grepout2 =  subprocess.check_output(cmdFls2, shell=True)
					except:
						True
						#print "not present in both 39 and V6_DBs"
						try:
							True
							cmdFls3 = "for i in /media/varma/My_Passport_1/Backup_Varma/RenduOsler2/work_WGS/V6_12345_*_variants_Varscan_SNPINDELS.vcf; do LANG=C grep -m 1 -w '"+str(infll[0]+'\t'+infll[1])+"' $i && break; done"
							grepout3 =  subprocess.check_output(cmdFls3, shell=True)
						except:
							#print infll[0]+":"+infll[1]+"-"+infll[2]
							file2.write(str(infl_nhr[count1]))
				count1 += 1
			file2.close()
			return None



# file1: Input positions SNP
# write1: Output file 

clF1 = SearchDB().Lsearch_UCSC_RG_LCR_HD_UM75(sys.argv[1],sys.argv[2])

SDB = SearchDB().Lsearch_WGS_WES(sys.argv[2],sys.argv[3])




