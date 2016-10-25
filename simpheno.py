#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from math import log10, sqrt, log, exp
from plinkio import plinkfile
from random import gauss, uniform, expovariate
from exceptions import *
from parser import *
from writer import *

class Simpheno():
	
	def __init__(self, inargs):
		#---------------------------------------------------------#
		self.description = "This shape has not been described yet"
		self.author = "Ilya Y. Zhbannikov"
		#---------------------------------------------------------#
		# Initialising class members #
		self.inType = inargs.itype
		self.oType = inargs.otype
		# Lists to save the results
		self.dtrait = None
		self.ctrait = None
		self.strait = None
		self.indiv_id = None
		# Paths #
		self.datapath = None
		self.cefile = None
		self.output_prefix = None
		# Genetics #
		self.genoMatrix = None
		self.alleleFreq = None
		self.snpeff = dict()
		# Misc #
		self.h = inargs.h
		self.alpha = inargs.alpha
		
	
	
	#----------------- Methods ------------------#
	
	def prepareCE(self):
		""" Prepares a list of effects from causal SNPs """
		
		eff = dict()
		
		# Read cefile
		f = open(self.cefile)
		lines = f.readlines()
		f.close()
		
		for l in lines :
			splits = l.split(':')
			eff[int(splits[0])] = float(splits[1].replace('\n',''))
		
		self.snpeff = eff

	
	def prepareMatrix( self, fname ):
		"""Prepares genotype matrix"""
		
		p = Parser()
		
		if self.inType == "plink" :
			self.genoMatrix = p.parse_plink(fname)
		elif self.inType == "ms" :
			self.genoMatrix = p.parse_ms(fname)
		elif self.inType == "genome" :
			self.genoMatrix = p.parse_genome(fname)
		
		""" Prepare genotype frequencies """
		freq = []
		for j in range(len(self.genoMatrix[0][0][0])) :
			z = 0.0
			o = 0.0
			t = 0.0
			for i in range(len(self.genoMatrix[0][0])) :
				snp = float(self.genoMatrix[0][0][i][j])
				if snp == 0 :
					z += 1.0
				elif snp == 1 :
					o += 1.0
				elif snp == 2 :
					t += 1.0
					
			nrow = len(self.genoMatrix[0][0])
			fA = z/nrow + 0.5*o/nrow
			fB = t/nrow + 0.5*o/nrow
			
			if fA == 0.0 :
				fA = 1e-12
				fB = 1 - fA
			if fB == 0.0 :
				fB = 1e-12
				fA = 1 - fB
		
			freq.append(min(fA,fB))
		
		self.alleleFreq = freq	
	
		
	def getId(self):
		"""Extracts ids from provided plink-formatted file"""
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			raise plinkioReadError
		
		sample_list = plink_file.get_samples( )
		if not sample_list :
			raise plinkioGetSamplesError
	
		ids = []
		for sample in sample_list:
			#print sample.fid, sample.iid, sample.phenotype
			ids.append([sample.fid, sample.iid])
	
		self.indiv_id = ids

	def simDichotomousPhe(self, cov, p0) :
		"""Simulates dichotomous phenotype"""
		
		bin_trait = []
		
		if cov != None:
			pass # TODO
		else :
			for row in self.genoMatrix[0][0] :
				sum_wij_ui = 0.0
				for g,freq,j in zip(row, self.alleleFreq, range(len(row))) : #, locus_list) :
					g = float(g)
					wij = 0.0
					if len(self.snpeff) == 0 :
						wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
						sum_wij_ui += g*wij
					else :
						if j in self.snpeff :
							sum_wij_ui += g*self.snpeff[j]
					#else :
					#	print "Warning: SNP with index", j, "not found."
			
				z = sum_wij_ui + gauss(0, 1)	
				pr = 1/(1+exp(-z))
			
				bt = 1 if pr > uniform(0,1) else 0
				
				#bin_trait.append([sample.fid, sample.iid, bt])
				bin_trait.append(bt)
	
		self.dtrait = bin_trait

	def simContinuousPhe(self, cov):
		"""Used to simulate continuous phenotype"""
		cont_trait = []
		
		if cov != None:
			pass # TODO
		else :
			for row in self.genoMatrix[0][0] :
				sum_wij_ui = 0.0
				for g,freq,j in zip(row, self.alleleFreq, range(len(row))) :
					g = float(g)
					wij = 0.0
					if len(self.snpeff) == 0 :
						wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
						sum_wij_ui += g*wij
					else :
						if j in self.snpeff :
							sum_wij_ui += g*self.snpeff[j]
			
				z = sum_wij_ui + gauss(0, 1)	
				
				cont_trait.append(z)
			
		self.ctrait = cont_trait

	
	def prepareCovariates(cov):
		covar = []
		mean_cov = prepareMeanCov(cov)
		for i in range(n):
			row = []
			for m in mean_cov :
				row.append(gauss(m, sqrt(m/10)))
			covar.append(row)
		return covar
	
	
	def simulWeib(self, cov, lambd, rho, rateC):
		"""Simulates survival time with Weibul distribution
			lamb = scale parameter in h0()
			rho = shape parameter in h0()
			rateC = rate parameter of the exponential distribution of C
		"""
		surv_trait = []
		
		for row in self.genoMatrix[0][0]:
			sum_wij_ui = 0.0
			for g,freq,j in zip(row,self.alleleFreq, range(len(row))) :
				g = float(g)
				wij = 0.0
				if len(self.snpeff) == 0 :
					wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					sum_wij_ui += g*wij
				else :
					if j in self.snpeff :
						sum_wij_ui += g*self.snpeff[j]
					
				
			v = uniform(0,1)
			
			Tlat = pow((-1*log(v) / (lambd * exp(sum_wij_ui))),(1/rho))
			# censoring times
			C = expovariate(rateC)
			# follow-up times and event indicators
  			time = min(Tlat, C)
  			status = 1 if Tlat >= C else 0
  			#surv_trait.append([time, status])
  			surv_trait.append(time)
  		self.strait = surv_trait
  	
  	
	def simulGomp(self, cov, lambd, rho, rateC):
		# lamb = scale parameter in h0()
		# rho = shape parameter in h0()
		# rateC = rate parameter of the exponential distribution of C
	
		# Gompertz latent event times
		surv_trait = []	
		
		for row in self.genoMatrix[0][0]:
			sum_wij_ui = 0.0
			for g,freq,j in zip(row,self.alleleFreq, range(len(row))) :
				g = float(g)
				wij = 0.0
				if len(self.snpeff) == 0 :
					wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					sum_wij_ui += g*wij
				else :
					if j in self.snpeff :
						sum_wij_ui += g*self.snpeff[j]
				
				
			v = uniform(0,1)
			
			Tlat = 1/self.alpha * log(1 - (self.alpha*log(v)) / (lambd * exp(sum_wij_ui)) )
			# censoring times
			C = expovariate(rateC)
			# follow-up times and event indicators
  			time = min(Tlat, C)
  			status = 1 if Tlat >= C else 0
  			surv_trait.append(time) #surv_trait.append([time, status])
  		self.strait = surv_trait
  		
	def prepare(self, path, cefile, outprefix) :
		self.datapath = path
		self.cefile = cefile
		self.output_prefix = outprefix
		
		# Prepare SNP effects from file
		if self.cefile != None :
			self.prepareCE()
	
		self.prepareMatrix(path)
		
	def saveData(self):
		
		if self.dtrait != None :
			# Saving dichotomous phenotype #
			ofname = self.output_prefix + "_pheno_bin.txt"
			f = open(ofname, "w")
			
			for d in self.dtrait :
				f.write(str(d))
				f.write('\n')
			f.close()
			
			# Formatted writing #
			self.writeInSpecificFormat(ofname, self.dtrait)
	
		if self.ctrait != None :
			# Saving continuous phenotype #
			ofname = self.output_prefix + "_pheno_cont.txt"
			f = open(ofname, "w")
			for c in self.ctrait :
				f.write(str(c))
				f.write('\n')
			f.close()
			
			# Formatted writing #
			self.writeInSpecificFormat(ofname, self.ctrait)
	
		if self.strait != None :
			# Saving continuous phenotype #
			ofname = self.output_prefix + "_pheno_surv.txt"
			f = open(ofname, "w")
			for s in self.strait :
				f.write(str(s))
				f.write('\n')
			f.close()
			
			# Formatted writing #
			self.writeInSpecificFormat(ofname, self.strait)
	
	
		if self.indiv_id != None :
			# Saving fid and id #
			ofname = self.output_prefix + "_id.txt"
			f = open(ofname, "w")
			for i in self.indiv_id :
				f.write('\t'.join(map(str,i)))
				f.write('\n')
			f.close()
			
	
		
		
	def writeInSpecificFormat(self, ofname, trait) :
		w = Writer()
		if self.oType == "emmax" :
			w.convert2emma(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
		elif self.oType == "plink" :
			w.convert2plink(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname, diploid=1)
		elif self.oType == "blossoc" :
			w.convert2blossoc(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
		elif self.oType == "qtdt" :
			w.convert2qtdt(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
		elif self.oType == "tassel" :
			w.convert2tassel(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
		
	