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

def simDichotomousPhe(path, freqpath, b0, cbetta, cov, p0) :
	# Read plink files:
	plink_file = plinkfile.open( path )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
	
	# Reading allele frequencies:
	f = open(freqpath, "rU")
	freq_lines = f.readlines()
	f.close()
	frequencies = []
	for line in freq_lines:
		if len(line.strip()) != 0:
			frequencies.append(float(line.split('\t')[1].split(' ')[0]))
	
	
	# Simulate binary (dichotomous) trait:
	bin_trait = []
	matrix = prepareMatrix(plink_file)
	if cov != None:
		for row, sample, c in zip(matrix, sample_list, cov ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
				for cc in c :
					csum = beta_p(cbetta, freq)*cc
				summary_genotype += csum
			z = float(b0) + summary_genotype + gauss(1, 0)
			#z = float(b0) + summary_genotype + sum(c) + gauss(1, 0)
			pr = 1/(1+exp(-z))
			bt = 1 if pr > p0 else 0
			#print pr, z, bt
			bin_trait.append([sample.fid, sample.iid, bt])
		#print len(plink_file)
		
	else :
		for row, sample in zip(matrix, sample_list ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
			z = float(b0) + summary_genotype + gauss(1, 0)
			pr = 1/(1+exp(-z))
			bt = 1 if pr > p0 else 0
			print pr, z, bt
			bin_trait.append([sample.fid, sample.iid, bt])
	
	return bin_trait

def simContinuousPhe(path, freqpath, b0, cbetta, cov):
	# Read plink files:
	plink_file = plinkfile.open( path )
	if not plink_file.one_locus_per_row( ):
		print( "This script requires that snps are rows and samples columns." )
		exit( 1 )
		
	sample_list = plink_file.get_samples( )
	locus_list = plink_file.get_loci( )
		
	#print "\n"
	# Reading allele frequencies:
	f = open(freqpath, "rU")
	freq_lines = f.readlines()
	f.close()
	frequencies = []
	for line in freq_lines:
		if len(line.strip()) != 0:
			frequencies.append(float(line.split('\t')[1].split(' ')[0]))
	
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
		for row, sample in zip( matrix, sample_list ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
			ct = float(b0) + summary_genotype + gauss(1, 0)
			cont_trait.append([sample.fid, sample.iid, ct])
	
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
	
def prepareCovariates(cov, n):
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
def simulWeib(lambd, rho, beta, rateC, covar, betas):
  # Weibull latent event times
  v = uniform(0,1)
  cb = 0
  for c, b in zip(covar, betas) :
  	 cb += c*b
  #print cb
  Tlat = pow((- log(v) / (lambd * exp(cb))),(1/rho))
  #print Tlat
  # censoring times
  C = expovariate(rateC)
  #print C
  # follow-up times and event indicators
  time = min(Tlat, C)
  status = 1 if Tlat <= C else 0
  return [time, status]
  
def simSurvPhe(lambd, rho, beta, rateC, covariates, betas) :
	result = []
	for c in covariates :
		result.append(simulWeib(lambd, rho, beta, rateC, c, betas))
	
	return result