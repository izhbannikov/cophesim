#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from math import log10, sqrt, log, exp
from plinkio import plinkfile
from random import gauss, uniform, expovariate
from exceptions import *
from parser import *

class Simpheno():
	
	h = 0.8
	alpha = 0.2138
	snpeff = dict()
	inType = "plink"
	
	# Lists to save the results
	dtrait = None
	ctrait = None
	strait = None
	indiv_id = None
	
	datapath = None
	cefile = None
	output_prefix = None
	
	genoMatrix = None
	
	def __init__(self, h=0.8, alpha=0.2138, intype="plink"):
		self.h = h
		self.alpha = alpha
		self.description = "This shape has not been described yet"
		self.author = "Ilya Y. Zhbannikov"
		print intype
		self.inType = intype
	
	#----------------- Methods ------------------#
	
	def prepareCE(self):
		""" Prepares a list of effects from causal SNPs """
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			exit( 1 )
		
		eff = dict()
		
		# Read cefile
		f = open(self.cefile)
		lines = f.readlines()
		f.close()
		
		for l in lines :
			splits = l.split(':')
			eff[splits[0]] = float(splits[1].replace('\n',''))
		
		self.snpeff = eff
		
		return eff

	
	def prepareMatrix( self, fname ):
		"""Prepares genotype matrix"""
		matrix = []
		p = Parser()
		print self.inType
		if self.inType == "plink" :
			matrix = p.parse_plink(fname)
		elif self.inType == "ms" :
			data = p.parse_ms(fname)
			matrix = data[0][0]
		elif self.inType == "genome" :
			data = p.parse_genome(fname)
			matrix = data[0][0]
		
		print matrix
		return matrix
	
	def getFreq( self, plink_file ):
		"""Calculates genotype frequency"""
		
		freq = []
		for col in plink_file:
			z = 0.0
			o = 0.0
			t = 0.0
			for i in range(len(col)) :
				if col[i] == 0 :
					z += 1.0
				elif col[i] == 1 :
					o += 1.0
				elif col[i] == 2 :
					t += 1.0
			fA = z/len(col) + 0.5*o/len(col)
			fB = t/len(col) + 0.5*o/len(col)
		
			if fA == 0.0 :
				fA = 1e-12
				fB = 1 - fA
			if fB == 0.0 :
				fB = 1e-12
				fA = 1 - fB
		
			freq.append(min(fA,fB)) 
	
		return freq
		
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
		
		# Read plink files:
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			raise plinkioReadError
		
		
		sample_list = plink_file.get_samples( )
		if not sample_list :
			raise getSamplesError
			
		locus_list = plink_file.get_loci( )
	
		# MAFs:
		frequencies = self.getFreq(plink_file)
		# Simulate binary (dichotomous) trait:
		bin_trait = []
		matrix = self.prepareMatrix(plink_file)
		
		
		if cov != None:
			pass # TODO
		else :
			for row, sample in zip(matrix, sample_list):
				sum_wij_ui = 0.0
				for g,freq,j, locus in zip(row,frequencies, range(len(row)), locus_list) :
					wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					
					if locus.name in self.snpeff :
						sum_wij_ui += wij*self.snpeff[locus.name]
			
				z = sum_wij_ui + gauss(0, 1)	
				pr = 1/(1+exp(-z))
			
				bt = 1 if pr > p0 else 0
				
				bin_trait.append([sample.fid, sample.iid, bt])
	
		self.dtrait = bin_trait

	def simContinuousPhe(self, cov):
		"""Used to simulate continuous phenotype"""
		
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			raise plinkioReadError
		
		sample_list = plink_file.get_samples( )
		if not sample_list :
			raise plinkioGetSamplesError
			
		locus_list = plink_file.get_loci( )
		if not locus_list :
			raise plinkioGetLocusError
			
		# MAFs:
		frequencies = self.getFreq(plink_file)
	
		# Simulate continuous trait:
		cont_trait = []
		matrix = self.prepareMatrix(plink_file)
		
		
		if cov != None:
			pass # TODO
		else :
			for row, sample in zip(matrix, sample_list):
				sum_wij_ui = 0.0
				for g,freq,j, locus in zip(row,frequencies, range(len(row)), locus_list) :
					wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					
					if locus.name in self.snpeff :
						sum_wij_ui += wij*self.snpeff[locus.name]
			
				z = sum_wij_ui + gauss(0, 1)	
				
				cont_trait.append([sample.fid, sample.iid, z])
			
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
		
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			raise plinkioReadError
		
		sample_list = plink_file.get_samples( )
		if not sample_list :
			raise plinkioGetSamplesError
			
		locus_list = plink_file.get_loci( )
		if not locus_list :
			raise plinkioGetLocusError
		
	
		# MAFs
		frequencies = self.getFreq(plink_file)
		# Prepare genotype matrix
		matrix = self.prepareMatrix(plink_file)
		

		# Weibull latent event times
		surv_trait = []
		
		for row, sample in zip(matrix, sample_list):
			sum_wij_ui = 0.0
			for g,freq,j, locus in zip(row,frequencies, range(len(row)), locus_list) :
				wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					
				if locus.name in self.snpeff :
					sum_wij_ui += wij*self.snpeff[locus.name]
				
			v = uniform(0,1)
			
			Tlat = pow((-1*log(v) / (lambd * exp(sum_wij_ui))),(1/rho))
			# censoring times
			C = expovariate(rateC)
			# follow-up times and event indicators
  			time = min(Tlat, C)
  			status = 1 if Tlat >= C else 0
  			surv_trait.append([sample.fid, sample.iid, time, status])
  		
  		self.strait = surv_trait
  	
  	
	def simulGomp(self, cov, lambd, rho, rateC):
		# lamb = scale parameter in h0()
		# rho = shape parameter in h0()
		# rateC = rate parameter of the exponential distribution of C
	
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			exit( 1 )
		
		sample_list = plink_file.get_samples( )
		locus_list = plink_file.get_loci( )
	
		# MAFs:
		frequencies = self.getFreq(plink_file)
		# Prepare genotype matrix
		matrix = self.prepareMatrix(plink_file)
	
		# Weibull latent event times
		surv_trait = []	
		
		for row, sample in zip(matrix, sample_list):
			sum_wij_ui = 0.0
			for g,freq,j, locus in zip(row,frequencies, range(len(row)), locus_list) :
				wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
					
				if locus.name in self.snpeff :
					sum_wij_ui += wij*self.snpeff[locus.name]
				
			v = uniform(0,1)
			
			Tlat = 1/self.halpha * log(1 - (self.halpha*log(v)) / (lambd * exp(sum_wij_ui)) )
			# censoring times
			C = expovariate(rateC)
			print Tlat, C
			# follow-up times and event indicators
  			time = min(Tlat, C)
  			status = 1 if Tlat >= C else 0
  			surv_trait.append([sample.fid, sample.iid, time, status])
  		
  		self.strait = surv_trait
  		
	def prepare(self, path, cefile, outprefix) :
		self.datapath = path
		self.cefile = cefile
		self.output_prefix = outprefix
		
		# Prepare SNP effects from file
		if self.cefile != None :
			self.prepareCE()
	
		# Get id of each individual
		#self.getId()
		
		self.prepareMatrix(path)
		
	def saveData(self):
	
		if self.dtrait != None :
			# Saving dichotomous phenotype:
			f = open(self.output_prefix + "_pheno_bin.txt", "w")
			for d in self.dtrait :
				f.write('\t'.join(map(str,d)))
				f.write('\n')
			f.close()
	
		if self.ctrait != None :
			# Saving continuous phenotype:
			f = open(self.output_prefix + "_pheno_cont.txt", "w")
			for c in self.ctrait :
				f.write('\t'.join(map(str,c)))
				f.write('\n')
			f.close()
	
		if self.strait != None :
			# Saving continuous phenotype:
			f = open(self.output_prefix + "_pheno_surv.txt", "w")
			for c in self.strait :
				f.write('\t'.join(map(str,c)))
				f.write('\n')
			f.close()
	
	
		if self.indiv_id != None :
			# Saving fid and id:
			f = open(self.output_prefix + "_id.txt", "w")
			for i in self.indiv_id :
				f.write('\t'.join(map(str,i)))
				f.write('\n')
			f.close()
	
			
	