#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from math import log10, sqrt, log, exp
from plinkio import plinkfile
from random import gauss, uniform, expovariate

def beta_p(cbetta, maf):
	res = cbetta*abs(log10(maf))
	return res
	
def prepareMatrix( plink_file ):
	matrix = []
	for row in plink_file:
		for i in range(len(row)):
			matrix.append([])
		break
	
	for row in plink_file:
		for i in range(len(row)):
			matrix[i].append(row[i])
	return matrix
	
def getFreq( plink_file ):
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

def simDichotomousPhe(idata, path, b0, cbetta, cov, p0) :
	# Read plink files:
	plink_file = plinkfile.open( idata )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	# MAFs:
	frequencies = getFreq(plink_file)
	# Simulate binary (dichotomous) trait:
	bin_trait = []
	matrix = prepareMatrix(plink_file)
	
	if cov != None:
		for row, sample, c in zip(matrix, sample_list, cov ):
			summary_genotype = 0
			#for g,freq in zip(row,frequencies) :
			#	summary_genotype += beta_p(cbetta, freq)*float(g)
			#	for cc in c :
			#		csum = beta_p(cbetta, freq)*cc
			#	summary_genotype += csum
			#summary_genotype = sqrt(1+pisum)*gauss(0,1)
			z = float(b0) + summary_genotype + gauss(1, 0)
			#z = float(b0) + summary_genotype + sum(c) + gauss(1, 0)
			pr = 1/(1+exp(-z))
			bt = 1 if pr > p0 else 0
			#print pr, z, bt
			bin_trait.append([sample.fid, sample.iid, bt])
		
	else :
		h = 0.8
		#ui = [0.000001, 1.5]
		ui = []
		[ui.append(0.000001) for i in range(len(matrix[0]))]
		
		for row, sample in zip(matrix, sample_list):
			sum_wij_ui = 0.0
			for g,freq,j in zip(row,frequencies, range(len(row))) :
				wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
				sum_wij_ui += wij*ui[j]
			
			#va = sum_wij_ui*(1.0/pow(h, 2) - 1.0)
			z = sum_wij_ui + gauss(0, 1)	
			pr = 1/(1+exp(-z))
			
			bt = 1 if pr > p0 else 0
			#print pr, z, bt
			
			bin_trait.append([sample.fid, sample.iid, bt])
	
	return bin_trait

def simContinuousPhe(idata, path, b0, cbetta, cov):
	plink_file = plinkfile.open( idata )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	# MAFs:
	frequencies = getFreq(plink_file)
	
	# Simulate continuous trait:
	cont_trait = []
	matrix = prepareMatrix(plink_file)
	if cov != None:
		for row, sample, c in zip( matrix, sample_list, cov ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
			ct = float(b0) + summary_genotype + sum(c) + gauss(1, 0)
			cont_trait.append([sample.fid, sample.iid, ct])
	else :
		h = 0.8
		#ui = [0.000001, 1.5]
		ui = []
		[ui.append(0.000001) for i in range(len(matrix[0]))]
		
		for row, sample in zip(matrix, sample_list):
			sum_wij_ui = 0.0
			for g,freq,j in zip(row,frequencies, range(len(row))) :
				wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
				sum_wij_ui += wij*ui[j]
			
			#va = sum_wij_ui*(1.0/pow(h, 2) - 1.0)
			z = sum_wij_ui + gauss(0, 1)	
			
			cont_trait.append([sample.fid, sample.iid, z])
			
	return cont_trait

def getDichotomousPhePlink(path):
	# Read plink files:
	plink_file = plinkfile.open( path )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	dichot_trait = []
	for sample in sample_list:
		#print sample.fid, sample.iid, sample.phenotype
		dichot_trait.append([sample.fid, sample.iid, sample.phenotype])
	
	return dichot_trait
	
def getId(path):
	# Read plink files:
	plink_file = plinkfile.open( path )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	
	ids = []
	for sample in sample_list:
		#print sample.fid, sample.iid, sample.phenotype
		ids.append([sample.fid, sample.iid])
	
	return ids
	
def prepareMeanCov(cov):
	mean_cov = cov.split(',')
	for i in range(len(mean_cov)):
		mean_cov[i] = float(mean_cov[i])
	return mean_cov
	
def prepareCovariates(cov):
	covar = []
	mean_cov = prepareMeanCov(cov)
	for i in range(n):
		row = []
		for m in mean_cov :
			row.append(gauss(m, sqrt(m/10)))
		covar.append(row)
	return covar
	
# lamb = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C
def simulWeib(idata, path, b0, cbetta, cov, lambd, rho, rateC):
	plink_file = plinkfile.open( idata )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	# MAFs:
	frequencies = getFreq(plink_file)
	# Prepare genotype matrix
	matrix = prepareMatrix(plink_file)
	
	# Weibull latent event times
	surv_trait = []
	#v = uniform(0,1)	
	ui = []
	[ui.append(0.00001) for i in range(len(matrix[0]))]	
	
	for row, sample in zip(matrix, sample_list):
		sum_wij_ui = 0.0
		for g,freq,j in zip(row,frequencies, range(len(row))) :
			wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
			sum_wij_ui += wij*ui[j]

		v = uniform(0,1)
		Tlat = pow((-1*log(v) / (lambd * exp(sum_wij_ui))),(1/rho))
		# censoring times
		C = expovariate(rateC)
		print Tlat, C
		# follow-up times and event indicators
  		time = min(Tlat, C)
  		status = 1 if Tlat >= C else 0
  		surv_trait.append([time, status])
  		
  	return surv_trait
  	
  	
# lamb = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C
def simulGomp(idata, path, b0, cbetta, cov, lambd, rho, rateC):
	plink_file = plinkfile.open( idata )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	# MAFs:
	frequencies = getFreq(plink_file)
	# Prepare genotype matrix
	matrix = prepareMatrix(plink_file)
	
	# Weibull latent event times
	surv_trait = []
	#v = uniform(0,1)	
	ui = []
	[ui.append(0.00001) for i in range(len(matrix[0]))]	
	
	for row, sample in zip(matrix, sample_list):
		sum_wij_ui = 0.0
		for g,freq,j in zip(row,frequencies, range(len(row))) :
			wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
			sum_wij_ui += wij*ui[j]

		v = uniform(0,1)
		alpha = 0.2138
		Tlat = 1/alpha * log(1 - (alpha*log(v)) / (lambd * exp(sum_wij_ui)) )
		# censoring times
		C = expovariate(rateC)
		print Tlat, C
		# follow-up times and event indicators
  		time = min(Tlat, C)
  		status = 1 if Tlat >= C else 0
  		surv_trait.append([time, status])
  		
  	return surv_trait