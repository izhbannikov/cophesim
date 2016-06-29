#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from itertools import chain

def saveData(ct, dt, ids, output_prefix, covariates):
	
	if ct != None :
		# Saving continuous phenotype:
		f = open(output_prefix + "_pheno_cont.txt", "w")
		for c in ct :
			f.write('\t'.join(map(str,c)))
			f.write('\n')
		f.close()
	
	if dt != None :
		# Saving dichotomous phenotype:
		f = open(output_prefix + "_pheno_bin.txt", "w")
		for d in dt :
			f.write('\t'.join(map(str,d)))
			f.write('\n')
		f.close()
	
	if ids != None :
		# Saving fid and id:
		f = open(output_prefix + "_id.txt", "w")
		for i in ids :
			f.write('\t'.join(map(str,i)))
			f.write('\n')
		f.close()
	
	# Saving covariates:
	if covariates != None :
		f = open(output_prefix + "_covariates.txt", "w")
		for i, c in zip(ids, covariates) :
			row = list(chain.from_iterable([i, c]))
			f.write('\t'.join(map(str,row)))
			f.write('\n')
			
		f.close()