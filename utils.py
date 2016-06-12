#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

def saveData(ct, dt, ids, output_prefix):
	# Saving continuous phenotype:
	f = open(output_prefix + "_pheno_cont.txt", "w")
	for c in ct :
		f.write('\t'.join(map(str,c)))
		f.write('\n')
	f.close()
	# Saving dichotomous phenotype:
	f = open(output_prefix + "_pheno_bin.txt", "w")
	for d in dt :
		f.write('\t'.join(map(str,d)))
		f.write('\n')
	f.close()
	# Saving fid and id:
	f = open(output_prefix + "_id.txt", "w")
	for i in ids :
		f.write('\t'.join(map(str,i)))
		f.write('\n')
	f.close()