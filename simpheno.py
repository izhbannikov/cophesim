#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from math import *
from plinkio import plinkfile
from random import gauss

def beta_p(cbetta, maf):
	res = cbetta*abs(log10(maf))
	return res


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
	if cov != None:
		for locus, row, sample, c in zip( locus_list, plink_file, sample_list, cov ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
			ct = float(b0) + summary_genotype + sum(c) + gauss(1, 0)
			cont_trait.append([sample.fid, sample.iid, ct])
	else :
		for locus, row, sample in zip( locus_list, plink_file, sample_list ):
			summary_genotype = 0
			for g,freq in zip(row,frequencies) :
				summary_genotype += beta_p(cbetta, freq)*float(g)
			ct = float(b0) + summary_genotype + gauss(1, 0)
			cont_trait.append([sample.fid, sample.iid, ct])
	
	return cont_trait

def simDichotomousPhe(path):
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