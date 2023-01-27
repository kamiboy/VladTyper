import os
import numpy as np
from time import time
from bedbug import BEDBUG

TEST_MODE = False

# Try to predict Sid Sda antigen type (untested)
PREDICT_SID = True

# Try to predict Lewis Lea and Leb antigen types (low  accuracy)
PREDICT_LEWIS = False

# directory containing bed/bim/fam genotypes
snpdir = '/data/preprocessed/genetics/dbds_freeze_20210503/'

# prefix of each bed/bim/fam chromosome file
chr_file = '20210503_chr'

# postfix of each bed/bim/fam chromosome file
chr_file_postfix = '_DBDS_FINAL'

# directory where all output files are written 
workdir = '/users/home/cammos/export/'

# directory where blood type data is read from
datadir = '/data/projects/rbc_gwas/'

# (optional) name of file with abo+rh serology data in compounded format e.g (AB+)
file_abo_rh = ''

# (optional) name of file (if any) with other antingen serology data, if any
file_antigens =''

# file describing how each antigen is predicted using a specific genotype
file_genotype_definitions = 'antigen.info.tsv'

def FBeta(B, precision, recall):
	return(((1 + (B**2)) * (precision * recall)) / ((B**2) * precision + recall))

class Variant:
	def __init__(self, record):
		self.chr = record[0]
		if not raw and self.chr == '23':
			self.chr = 'X'
		self.pos = int(record[1])
		self.id = record[2]
		self.ref = record[3]
		self.alt = record[4]

class Genotypes:
	def __init__(self, antigens):
		self.antigens = antigens 
		self.variants = {}

	def type(self, id, type):
		if id in self.variants:
			return
		self.variants[id] = type

class VladGenotyper:
	def __init__(self, phenotypes):
		self.donors = {}
		self.phenotypes = phenotypes

	#--- Lea/Leb
	def predictLewis(self):
		print("Predicting Lewis using genotypes")
		if not TEST_MODE:
			report = open(workdir + 'predictions.report.csv', 'a')

		chr = '19'
		if chr_file_postfix == None:
			bed = BEDBUG(snpdir+chr_file)
		else:
			bed = BEDBUG(snpdir+chr_file+chr+chr_file_postfix)
		
		phenotyped = 0
		congruent = 0
		genotyped = 0
		TP = 0
		TN = 0
		FP = 0
		FN = 0
		infoscore = 'nan'
		freq = 'nan'

		variant1 = None
		(variants, cases, _, rs601338) = bed.region(chr, 48703417, 48703417, [])
		if len(variants) > 0:
			variant1 = variants[0]
			if variant1.allele1 == 'G' and variant1.allele2 == 'A':
				variant1.allele1 = 'A'
				variant1.allele2 = 'G'
				for index in range(0,len(cases)):
					if rs601338[index] == 0:
						rs601338[index] = 2
					elif rs601338[index] == 2:
						rs601338[index] = 0

		if variant1 == None or not (variant1.allele1 == 'A' and variant1.allele2 == 'G'):
			variant1 = None
			rs601338 = [-1]*len(cases)

		variant2 = None
		(variants, _, _, rs3745635) = bed.region(chr, 5844332, 5844332, [])
		if len(variants) > 0:
			variant2 = variants[0]
			if variant2.allele1 == 'C' and variant2.allele2 == 'T':
				variant2.allele1 = 'T'
				variant2.allele2 = 'C'
				for index in range(0,len(cases)):
					if rs3745635[index] == 0:
						rs3745635[index] = 2
					elif rs3745635[index] == 2:
						rs3745635[index] = 0

		if variant2 == None or not (variant2.allele1 == 'T' and variant2.allele2 == 'C'):
			variant2 = None
			rs3745635 = [-1]*len(cases)

		variant3 = None
		(variants, _, _, rs28362459) = bed.region(chr, 5844781, 5844781, [])
		if len(variants) > 0:
			variant3 = variants[0]
			if variant3.allele1 == 'A' and variant3.allele2 == 'C':
				variant3.allele1 = 'C'
				variant3.allele2 = 'A'
				for index in range(0,len(cases)):
					if rs28362459[index] == 0:
						rs28362459[index] = 2
					elif rs28362459[index] == 2:
						rs28362459[index] = 0

		if variant3 == None or not (variant3.allele1 == 'C' and variant3.allele2 == 'A'):
			variant3 = None
			rs28362459 = [-1]*len(cases)

		if variant1 == None or variant2 == None or variant3 == None:
			print('Warning: Lewis prediction not possible due to missing genotypes')

		leapos = 0
		lebpos = 0
		leneg = 0
		ambiguous = 0
		unknown = 0
		for index in range(0,len(cases)):
			id = cases[index]

			if rs3745635[index] == 0 or rs28362459[index] == 0:
				FUT3 = False
			elif (rs3745635[index] == 2 and rs28362459[index] == 2) or (rs3745635[index] == 1 and rs28362459[index] == 2) or (rs3745635[index] == 2 and rs28362459[index] == 1):
				FUT3 = True
			else:
				FUT3 = None 

			if rs601338[index] == 2 or rs601338[index] == 1:
				FUT2 = True
			elif rs601338[index] == 0:
				FUT2 = False
			else:
				FUT2 = None

			if FUT2 == None:
				unknown += 1
			elif (rs3745635[index] == 1 and rs28362459[index] == 1):
				ambiguous += 1
			elif FUT3 == None:
				unknown += 1
			elif FUT3 != None and FUT2 != None:
				if id not in self.donors:
					donor = Donor(id)
					self.donors[id] = donor

				donor = self.donors[id]
			
			phenotype = None
			if FUT3 == False and (FUT2 == False or FUT2 == True):
				leneg += 1
				donor.type('Lea', False, 'rs7224888/rs3745635/rs601338')
				donor.type('Leb', False, 'rs7224888/rs3745635/rs601338')
				phenotype = {'Lea':False,'Leb':False}
				#phenotype = {'Lea':False}
				#phenotype = {'Leb':False}
			elif FUT3 == True and FUT2 == False:
				leapos += 1
				donor.type('Lea', True, 'rs7224888/rs3745635/rs601338')
				donor.type('Leb', False, 'rs7224888/rs3745635/rs601338')
				phenotype = {'Lea':True,'Leb':False}
				#phenotype = {'Lea':True}
				#phenotype = {'Leb':False}
			elif FUT3 == True and FUT2 == True:
				lebpos += 1
				donor.type('Lea', False, 'rs7224888/rs3745635/rs601338')
				donor.type('Leb', True, 'rs7224888/rs3745635/rs601338')
				phenotype = {'Lea':False,'Leb':True}
				#phenotype = {'Lea':False}
				#phenotype = {'Leb':True}

			if phenotype != None:
				genotyped += 1
				if donor.id in self.phenotypes.donors:
					serologi = self.phenotypes.donors[donor.id]
					for antigen in phenotype:
						if antigen in serologi.bloodtypes:
							phenotyped = phenotyped + 1
							a = serologi.bloodtypes[antigen]
							b = donor.bloodtypes[antigen]

							if a > 0:
								if b > 0:
									TP = TP + 1
								if b < 0:
									FN = FN + 1
									donor.conflicts.add(antigen)
							if a < 0:
								if b < 0:
									TN = TN + 1
								if b > 0:
									donor.conflicts.add(antigen)
									FP = FP + 1

							if np.sign(a) == np.sign(b):
								congruent = congruent + 1

		if genotyped > 0 and phenotyped > 0:
			acc = str(round(float(congruent*100)/float(phenotyped),6))+'%'
		else:
			acc = 'NA'

		precision = float('nan')
		sensitivity = float('nan')
		specificity = float('nan')

		if TP > 0 or FP > 0:
			precision = float(TP) / float(TP+FP)
		if TP > 0 or FN > 0:
			sensitivity = float(TP) / float(TP+FN)
		if TN > 0 or FP > 0:
			specificity = float(TN) / float(TN+FP)

		fscore = 'NA'
		if precision > 0:
			fscore = str(FBeta(1.0, precision, sensitivity))

		print('Lewis: Le(a+b-): %i (%.2f), Le(a-b+): %i (%.2f), Le(a-b-): %i (%.2f), ambiguous: %i (%.2f), unknown: %i (%.2f), acc: %i/%i (%s), TPR: %.4f, TNR: %.4f, Fscore: (TP|%i TN|%i FP|%i FN|%i)=%s'%(leapos,leapos/(leapos+lebpos+leneg),lebpos,lebpos/(leapos+lebpos+leneg),leneg,leneg/(leapos+lebpos+leneg),ambiguous,ambiguous/len(cases),unknown,unknown/len(cases),congruent,phenotyped,acc,sensitivity,specificity,TP,TN,FP,FN,fscore ))
		if not TEST_MODE:
			report.write('Lewis\tLea/Leb\tNA\t'+str(genotyped)+'\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t'+str(phenotyped)+'\t'+str(congruent)+'\t'+acc+'\t'+str(sensitivity)+'\t'+str(specificity)+'\t'+str(TP)+'\t'+str(TN)+'\t'+str(FP)+'\t'+str(FN)+'\t'+fscore+'\n')
			report.close()

	#--- Sda
	def predictSid(self):
		print("Predicting Sda using genotypes")
		if not TEST_MODE:
			report = open(workdir + 'predictions.report.csv', 'a')

		chr = '17'

		if chr_file_postfix == None:
			bed = BEDBUG(snpdir+chr_file)
		else:
			bed = BEDBUG(snpdir+chr_file+chr+chr_file_postfix)
		
		phenotyped = 0
		congruent = 0
		genotyped = 0
		TP = 0
		TN = 0
		FP = 0
		FN = 0
		infoscore = 'nan'
		freq = 'nan'

		variant1 = None
		(variants, cases, _, rs7224888) = bed.region(chr, 49168801, 49168801, [])
		if len(variants) > 0:
			variant1 = variants[0]
			if variant1.allele1 == 'C' and variant1.allele2 == 'T':
				variant1.allele1 = 'T'
				variant1.allele2 = 'C'
				for index in range(0,len(cases)):
					if rs7224888[index] == 0:
						rs7224888[index] = 2
					elif rs7224888[index] == 2:
						rs7224888[index] = 0

		if variant1 == None or not (variant1.allele1 == 'T' and variant1.allele2 == 'C'):
			variant1 = None
			rs7224888 = [-1]*len(cases)

		variant2 = None
		(variants, _, _, rs72835417) = bed.region(chr, 49164280, 49164280, [])
		if len(variants) > 0:
			variant2 = variants[0]
			if variant2.allele1 == 'A' and variant2.allele2 == 'G':
				variant2.allele1 = 'G'
				variant2.allele2 = 'A'
				for index in range(0,len(cases)):
					if rs72835417[index] == 0:
						rs72835417[index] = 2
					elif rs72835417[index] == 2:
						rs72835417[index] = 0

		if variant2 == None or not (variant2.allele1 == 'G' and variant2.allele2 == 'A'):
			variant2 = None
			rs72835417 = [-1]*len(cases)

		if variant1 == None or variant2 == None:
			print('Warning: Sda prediction not possible due to missing genotypes')

		pos = 0
		neg = 0
		ambiguous = 0
		unknown = 0
		for index in range(len(cases)):
			id = cases[index]
			genotype = None

			if rs7224888[index] == 2 or rs72835417[index] == 2:
				genotype = False
				neg += 1
			elif (rs7224888[index] == 0 and rs72835417[index] == 0) or (rs7224888[index] == 0 and rs72835417[index] == 1) or (rs7224888[index] == 1 and rs72835417[index] == 0):
				genotype = True
				pos += 1
			elif rs7224888[index] == 1 and rs72835417[index] == 1:
				ambiguous += 1
			else:
				unknown += 1

			if genotype != None:
				genotyped += 1
				if id not in self.donors:
					donor = Donor(id)
					self.donors[id] = donor

				donor = self.donors[id]
				donor.type('Sda', genotype, 'rs7224888/rs72835417')

				if donor.id in self.phenotypes.donors:
					serologi = self.phenotypes.donors[donor.id]
					if 'Sda' in serologi.bloodtypes:
						phenotyped += 1
						a = serologi.bloodtypes['Sda']
						b = donor.bloodtypes['Sda']

						if a > 0:
							if b > 0:
								TP = TP + 1
							if b < 0:
								FN = FN + 1
								donor.conflicts.add('Sda')
						if a < 0:
							if b < 0:
								TN = TN + 1
							if b > 0:
								FP = FP + 1
								donor.conflicts.add('Sda')

						if np.sign(a) == np.sign(b):
							congruent = congruent + 1

		if genotyped > 0 and phenotyped > 0:
			acc = str(round(float(congruent*100)/float(phenotyped),6))+'%'
		else:
			acc = 'NA'

		precision = float('nan')
		sensitivity = float('nan')
		specificity = float('nan')

		if TP > 0 or FP > 0:
			precision = float(TP) / float(TP+FP)
		if TP > 0 or FN > 0:
			sensitivity = float(TP) / float(TP+FN)
		if TN > 0 or FP > 0:
			specificity = float(TN) / float(TN+FP)

		fscore = 'NA'
		if precision > 0:
			fscore = str(FBeta(1.0, precision, sensitivity))

		print('Sid Sda: pos: %i (%.2f), neg: %i (%.2f), ambiguous: %i (%.2f), unknown: %i (%.2f), acc: %i/%i (%s), TPR: %.4f, TNR: %.4f, Fscore: (TP|%i TN|%i FP|%i FN|%i)=%s'%(pos,pos/(pos+neg),neg,neg/(pos+neg),ambiguous,ambiguous/len(cases),unknown,unknown/len(cases),congruent,phenotyped,acc,sensitivity,specificity,TP,TN,FP,FN,fscore ))
		if not TEST_MODE:
			report.write('Sid\tSda\tNA\t'+str(genotyped)+'\tNA\tNA\tNA\tNA\tNA\tNA\t'+str(pos)+'\tNA\t'+str(phenotyped)+'\t'+str(congruent)+'\t'+acc+'\t'+str(sensitivity)+'\t'+str(specificity)+'\t'+str(TP)+'\t'+str(TN)+'\t'+str(FP)+'\t'+str(FN)+'\t'+fscore+'\n')
			report.close()

	def predictAntigens(self, file):
		print("Predicting phenotypes using genotypes")

		if not TEST_MODE:
			report = open(workdir + 'predictions.report.csv', 'wt')
			report.write('system\tantigens\tid\tgenotyped\tmissing\tinfoscore\tzygosity.ref\tzygosity.het\tzygosity.alt\tmaf\tpos.ref\tpos.alt\tphenotyped\tcongruent\tacc\tsensitivity\tspecificity\tTP\tTN\tFP\tFN\tFscore\n')

		if chr_file_postfix == None:
			bedfile = BEDBUG(snpdir+chr_file)

		file = open(file,'r')
		for line in file:
			if line[0] == '!':
				break
			if line == '\n' or line[0] == '#':
				continue
			items = line.strip('\n').split('\t')
			chr = items[0]
			pos = items[1]
			snp = items[3]
			ref = items[4]
			alt = items[5]
			antigen1 = items[6]
			antigen2 = items[7]
			system = items[8]

			reverse = False
			if chr_file_postfix != None:
				bedfile = BEDBUG(snpdir+chr_file+chr[3:]+chr_file_postfix)

			(variant, ids, genotypes) = bedfile.variant(chr[3:], int(pos), ref, alt, [])

			if variant == None or ids == None or genotypes == None:
				(variant, ids, genotypes) = bedfile.variant(chr[3:], int(pos), alt, ref, [])
				reverse = True

			if variant == None or ids == None or genotypes == None:
				print('Warning: ' + system + ' antigens ' + antigen1 + '/' + antigen2 + ' (' + snp + ') not found at ' + chr + ':' + pos + ':' + ref + '/' + alt)
			else:
				phenotyped = 0
				congruent = 0
				genotyped = 0
				positive1 = 0
				positive2 = 0
				TP = 0
				TN = 0
				FP = 0
				FN = 0

				infoscore = 'nan'
				freq = 'nan'

				for i in range(len(ids)):
					id = ids[i]
					genotype = genotypes[i]

					if reverse:
						if genotype != -1 and genotype != 1:
							if genotype == 0:
								genotype = 2
							elif genotype == 2:
								genotype = 0
							else:
								print('Error: Unknown genotype ' + str(genotype))
								exit(1)

					if genotype != -1:
						genotyped = genotyped + 1
						if id not in self.donors:
							donor = Donor(id)
							self.donors[id] = donor

						donor = self.donors[id]

						phenotype = None

						if antigen1 == '!' or antigen2 == '!':
							if (not reverse and genotype == 0) or (reverse and genotype == 2):
								if antigen1 != '!':
									donor.type(antigen1, False, snp)
									phenotype = {antigen1:False}
								elif antigen2 != '!':
									donor.type(antigen2, False, snp)
									phenotype = {antigen2:False}
								else:
									print('Error: Unknown genotype ' + str(genotype))
									exit(1)
							else:
								if antigen1 != '!':
									donor.type(antigen1, True, snp)
									phenotype = {antigen1:True}
									positive1 = positive1 + 1
								elif antigen2 != '!':
									donor.type(antigen2, True, snp)
									phenotype = {antigen2:True}
									positive2 = positive2 + 1
								else:
									print('Error: Unexpected antigens ' + antigen1 + '/' + antigen2)
									exit(1)
						elif antigen1 == '#' or antigen2 == '#':
							if (not reverse and genotype == 2) or (reverse and genotype == 0):
								if antigen1 != '#':
									donor.type(antigen1, True, snp)
									phenotype = {antigen1:True}
									positive1 = positive1 + 1
								elif antigen2 != '#':
									donor.type(antigen2, True, snp)
									phenotype = {antigen2:True}
									positive2 = positive2 + 1
								else:
									print('Error: Unknown genotype ' + str(genotype))
									exit(1)
							else:
								if antigen1 != '#':
									donor.type(antigen1, False, snp)
									phenotype = {antigen1:False}
								elif antigen2 != '#':
									donor.type(antigen2, False, snp)
									phenotype = {antigen2:False}
								else:
									print('Error: Unexpected antigens ' + antigen1 + '/' + antigen2)
									exit(1)
						else:
							if genotype == 0:
								donor.type(antigen1, True, snp)
								donor.type(antigen2, False, snp)
								phenotype = {antigen1:True,antigen2:False}
								positive1 = positive1 + 1
							elif genotype == 2:
								donor.type(antigen1, False, snp)
								donor.type(antigen2, True, snp)
								phenotype = {antigen1:False,antigen2:True}
								positive2 = positive2 + 1
							elif genotype == 1:
								donor.type(antigen1, True, snp)
								donor.type(antigen2, True, snp)
								phenotype = {antigen1:True,antigen2:True}
								positive1 = positive1 + 1
								positive2 = positive2 + 1
							else:
								print('Error: Unexpected genotype ' + str(genotype))
								exit(1)

						if donor.id in self.phenotypes.donors:
							serologi = self.phenotypes.donors[donor.id]
							for antigen in phenotype:
								if antigen in serologi.bloodtypes:
									phenotyped = phenotyped + 1
									a = serologi.bloodtypes[antigen]
									b = donor.bloodtypes[antigen]
									
									if a > 0:
										if b > 0:
											TP = TP + 1
										if b < 0:
											FN = FN + 1
											donor.conflicts.add(antigen)
									if a < 0:
										if b < 0:
											TN = TN + 1
										if b > 0:
											FP = FP + 1
											donor.conflicts.add(antigen)

									if np.sign(a) == np.sign(b):
										congruent = congruent + 1

				if phenotyped > 0:
					acc = str(round(float(congruent*100)/float(phenotyped),6))+'%'
				else:
					acc = 'NA'

				precision = float('nan')
				sensitivity = float('nan')
				specificity = float('nan')

				if TP > 0 or FP > 0:
					precision = float(TP) / float(TP+FP)
				if TP > 0 or FN > 0:
					sensitivity = float(TP) / float(TP+FN)
				if TN > 0 or FP > 0:
					specificity = float(TN) / float(TN+FP)

				fscore = 'NA'
				if precision > 0:
					fscore = str(FBeta(1.0, precision, sensitivity))
				
				print('System ' + system + ' antigens ' + antigen1 + '/' + antigen2 + ' ' + snp + ' genotyped: ' + str(genotyped) + ', missing: ' + str(variant.missing) + ', zygosity: (ref|'+ str(variant.homREF) + ' het|' + str(variant.het) + ' alt|' + str(variant.homALT) + '), infoscore: ' + infoscore + ', freq: ' + freq + ', maf: ' + str(variant.maf) + ', +: ' + str(round(float(positive1*100)/float(genotyped),4)) +'%/' + str(round(float(positive2*100)/float(genotyped),4)) + '%, acc: ' + str(congruent) +'/' + str(phenotyped) + ' (' + acc + ')' + ', TPR: ' + str(sensitivity) + ', TNR: ' + str(specificity) + ', Fscore: (TP|'+str(TP)+' TN|'+str(TN)+' FP|'+str(FP)+' FN|'+str(FN)+')=' + fscore ) #', orientation: ' + orientation +
				if not TEST_MODE:
					report.write(system + '\t' + antigen1 + '/' + antigen2 + '\t' + snp + '\t' + str(genotyped) + '\t' + str(variant.missing) + '\t' + infoscore + '\t'+ str(variant.homREF) + '\t' + str(variant.het) + '\t' + str(variant.homALT) + '\t' + str(variant.maf) + '\t' + str(float(positive1)/float(genotyped)) +'\t' + str(float(positive2*100)/float(genotyped)) + '\t' + str(phenotyped) + '\t' + str(congruent) +'\t' + acc + '\t' + str(sensitivity) + '\t' + str(specificity) + '\t' + str(TP) + '\t' + str(TN) + '\t' + str(FP) + '\t' + str(FN) + '\t' + fscore + '\n')

		file.close()
		if not TEST_MODE:
			report.close()

	def exportConflicts(self, conflicts,badsamples):
		badsamples = open(badsamples, "wt")
		conflicts = open(conflicts, "wt")
		conflicts.write('id\tconflicts\tantigens\n')
		badsamples.write('id\n')

		mismatch = []
		for id in self.donors:
			donor = self.donors[id]
			if len(donor.conflicts) > 4:
				badsamples.write(id + '\n')

			if len(donor.conflicts) > 0:
				conflicts.write(donor.id + '\t' + str(len(donor.conflicts)) + '\t')
				first = True
				for antigen in donor.conflicts:
					if not first:
						conflicts.write('/')
					else:
						first = False
					conflicts.write(antigen)
				conflicts.write('\n')
		conflicts.close()
		badsamples.close()

class Donor:
	def __init__(self, id):
		self.id = id
		self.bloodtypes = {}
		self.sources = {}
		self.conflicts = set()
		self.typed = {}
		self.genotypes = {}

	def genotype(self, antigens, variant, type):

		if antigens not in self.genotypes:
			self.genotypes[antigens] = Genotypes(antigens)

		genotype = self.genotypes[antigens]
		id = variant.chr+':'+str(variant.pos)+':'+variant.ref+'/'+variant.alt
		genotype.type(id, type)
		id = variant.chr+':'+str(variant.pos)+':'+variant.alt+'/'+variant.ref
		if type == 2:
			type = 0
		elif type == 0:
			type = 2
		genotype.type(id, type)

	def type(self, antigen, status, source):
		conflict = False
		mismatch = False
		if antigen == '*':
			return conflict

		if status == False:
			status = -1
		elif status == True:
			status = 1
		else:
			print("Error: invalid blood type status '" + str(status) + "'")
			exit(1)

		if antigen in self.bloodtypes:
			self.typed[antigen] = self.typed[antigen] + 1
			self.sources[antigen].add(source)

			if (self.bloodtypes[antigen]< 0 and status > 0) or (self.bloodtypes[antigen] > 0 and status < 0):
				conflict = True
			self.bloodtypes[antigen] = self.bloodtypes[antigen] + status

		else:
			self.bloodtypes[antigen] = status
			self.typed[antigen] = 1
			self.sources[antigen] = {source}

		return conflict

class VladPhenotyper:
	def __init__(self):
		self.donors = {}
		self.phenotypes = 0
		self.antigens = {'A','B','D'}

	def load_ABO_Rh_types(self, file):
		file = open(file, "r")
		linenr = 0
		for line in file:
			linenr += 1

			if linenr == 1:
				continue

			items = line.rstrip().split('\t')

			id = items[0]
			type = items[-1]

			A = False
			B = False
			D = False

			if 'A' in type:
				A = True
			if 'B' in type:
				B = True
			if '+' in type:
				D = True

			if id in self.donors:
				donor = self.donors[id]
			else:
				donor = Donor(id)
				self.donors[id] = donor

			donor.type('A', A, "serologi")
			donor.type('B', B, "serologi")
			donor.type('D', D, "serologi")
		self.phenotypes = self.phenotypes + linenr - 1
		print("Number of phenotyped cases: " + str(len(self.donors)))
		print("Number of phenotypes: " + str(linenr-1))

	def load_antigen_types(self, file):
		file = open(file, "r")
		linenr = 0
		for line in file:
			linenr += 1

			if linenr == 1:
				continue

			items = line.rstrip().split('\t')

			id = items[0]
			antigen = items[1]
			if items[2] == '1':
				status = True
			elif items[2] == '0':
				status = False
			else:
				print('Error: Unknown phenotype status %s'%items[2])
				exit(1)

			self.antigens.add(antigen)

			if id in self.donors:
				donor = self.donors[id]
			else:
				donor = Donor(id)
				self.donors[id] = donor

			donor.type(antigen, status, "serologi")
		self.phenotypes = self.phenotypes + linenr - 1
		print("Number of phenotyped cases: " + str(len(self.donors)))
		print("Number of phenotypes: " + str(linenr-1))

def dump(serologi, genotype, file):
	print("Taking a dump:")
	ids = set()
	types = set()

	for id in serologi.donors:
		ids.add(id)
		donor = serologi.donors[id]
		for type in donor.bloodtypes:
			types.add(type)

	for id in genotype.donors:
		ids.add(id)
		donor = genotype.donors[id]
		for type in donor.bloodtypes:
			types.add(type)

	file = open(file, "wt")
	file.write('id\tantigen\tserologi\tgenotype\n')

	for type in types:
		for id in ids:
			s = 'NA'
			g = 'NA'

			if id in serologi.donors:
				donor = serologi.donors[id]
				if type in donor.bloodtypes:
					s = str(int(donor.bloodtypes[type]))

			if id in genotype.donors:
				donor = genotype.donors[id]
				if type in donor.bloodtypes:
					g = str(int(donor.bloodtypes[type]))

			if s != 'NA' or g != 'NA':
				file.write(id + '\t' + type + '\t' + s + '\t' + g + '\n')
	file.close()

timestamp = time()
phenotyper = VladPhenotyper()

if file_abo_rh != '':
	phenotyper.load_ABO_Rh_types(datadir+file_abo_rh)

if file_antigens != '':
	phenotyper.load_antigen_types(datadir+file_antigens)

print("Process took %0.3fs" % (time() - timestamp))

timestamp = time()
genotyper = VladGenotyper(phenotyper)

genotyper.predictAntigens(datadir+file_genotype_definitions)

if PREDICT_SID:
	genotyper.predictSid()

if PREDICT_LEWIS:
	genotyper.predictLewis()

genotyper.exportConflicts(workdir+'predictions.conflicts.tsv',workdir+'predictions.badsamples.tsv')

print("Process took %0.3fs" % (time() - timestamp))

if not TEST_MODE:
	timestamp = time()
	dump(phenotyper, genotyper, workdir+'predictions.tsv')
	print("Process took %0.3fs" % (time() - timestamp))

print('Done!')