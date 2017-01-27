#!/usr/bin/env python


import sys
import subprocess


#############------------

"""
Script to add gene expression FPKM scores from Gene-Epression-BodyMap2.0 public resource for 16 different tissues, annotation for each gene can be performed


INPUT and OUTPUT file for the script:

python Gene-Epression-BodyMap2.0.py /media/varma/My_Passport_2/Scripts-for-LocalDBsearch/OUT_DATA/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_GATK_ANNOVAR_het_PASS-1%_notinUCSC-inRG_notinLocalDB_noINDELs_inBAMcheck-exc-chr5-6_bed12-allScore2_UCSC.bed /media/varma/My_Passport_2/Gene-Expression_GETx/OUT_DATA/HCM_B00H7EW-B00H7EX-B00H7EY_chr1-MT_Variants_GATK_ANNOVAR_het_PASS-1%_notinUCSC-inRG_notinLocalDB_noINDELs_inBAMcheck-exc-chr5-6_bed12-allScore2_UCSC_GeneExp.bed

#Gene_Names	adipose	adrenal	blood	brain	breast	colon	heart	kidney	liver	lung	lymph	ovary	prostate	skeletal_muscle	testes	thyroid

"""

#####----------------

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
		self.ident_ranges_HMBM = defaultdict(list)

	def Search_ExpDB(self, GName):
		self.list1 = []
		DIR="/media/varma/My_Passport_2/Gene-Expression_GETx/FPKMs_Human-BodyMap2.0_mRNA-seq-data_gene.matrix.csv"
		try:
			True
			cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(DIR)+""
			grepout1 =  subprocess.check_output(cmdFls1, shell=True)
			self.list1 = grepout1.strip().split(",")
		except:
			True
			#print GName+" not present in DB"
			self.list1 = [str('"'+GName+'"'), '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
		return self.list1
	
	def funczip(self,lt1,lt2):
		return '\t'.join(i + ', ' + j for i,j in zip(lt1,lt2))

	def Human_BodyMap20(self,readfl1,writefl3):
		"""
		Calling HBM localsearch method
		"""
		with open(readfl1,'r') as f1, open(writefl3,'w') as output:
			first_line = f1.readline().strip()
			HeadDB = "Gene_Names	Adipose	Adrenal	Blood	Brain	Breast	Colon	Heart	Kidney	Liver	Lung	Lymph	Ovary	Prostate	Skeletal_muscle	Testes	Thyroid"
			output.write(first_line+"\t"+HeadDB+"\n")
			for line in f1:
				lineS = line.strip()
				value = line.strip().split()[33]
				val1 = value.split(",")
				if len(val1) == 1:
					ExpG = "\t".join(self.Search_ExpDB(val1[0]))
					output.write(lineS+"\t"+ExpG+"\n")
				elif len(val1) != 1:
					ExpG1 = self.Search_ExpDB(val1[0])
					ExpG2 = self.Search_ExpDB(val1[1])
					MergExpG = self.funczip(ExpG1, ExpG2)
					output.write(lineS+"\t"+MergExpG+"\n")
		print "Done seach for Exp-data of Human_BodyMap20...."
		return None



# file1: Input positions SNP
# write1: Output file 

clF1 = SearchDB().Human_BodyMap20(sys.argv[1],sys.argv[2])






