#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

import subprocess
import sys
import argparse
import os
import ntpath

# Importing custom modules
from simpheno import *
from utils import *

# Starts the simulation
def main() :

	#------------Setting-up options for argument parser -------------#
	usage = "cophesim.py -i <path to plink files> -o <output prefix> [other options]"
	parser = argparse.ArgumentParser(usage=usage)
	
	# Input parameters
	parser.add_argument("-i", "--input", action="store", type=str, dest="idata", default=None, required=True, help="Path with prefix to your input file(s).")
	parser.add_argument("-o", "--output", action="store", type=str, dest="output_prefix", default=None, required=True, help="Output prefix.")
	parser.add_argument("-type", action="store", type=str, dest="intype", default="plink", required=False, help="Input type: plink (for Plink, default), ms (for ms, msms, msHot), genome (for Genome).")
	
	
	# Trait (outcome) type
	parser.add_argument("-d", "--dichotomous",action="store_true", dest="dflag", default=True, required=False, help="A flag for dichotomous phenotype, True by default.")
	parser.add_argument("-c", "--continuous",action="store_true", dest="cflag", default=False, required=False, help="A flag for continous phenotype, False by default.")
	parser.add_argument("-s", "--suvival", action="store_true", dest="sflag", default=False, required=False, help="A flag to simulate survival phenotype, False by default.")
	
	# Effects from causal SNPs
	parser.add_argument("-ce", action="store", dest="ceff", default=None, required=False, help="A path to the file with effect of each causal SNP. Must be in format: snp_name:effect.")
	
	# Other parameters
	parser.add_argument("-cov", "--covariates", action="store", type=str, dest="cov", default=None, required=False, help="Mean values of covariates, must be enumerated with comma and no spaces")
	parser.add_argument("-p0", "--p0",action="store", type=float, dest="p0", default=0.5, required=False, help="Probability for logistic model.")
	parser.add_argument("-weib", action="store_true", dest="weib", default=True, required=False, help="A flag to use Weibull distribution for survival phenotype. True by default.")
	parser.add_argument("-gomp", action="store_true", dest="gomp", default=False, required=False, help="A flag to use Gompertz distribution for survival phenotype. False by default.")
	
	
	#--------- Checking the input arguments ----------#
	args = parser.parse_args()
	
	if not os.path.exists(os.path.dirname(args.idata)):
		print "Directory", args.indata, "not exists. Aborting."
		sys.exit(-1)
	
	if not os.path.exists(os.path.dirname(args.output_prefix)):
		print "Directory", args.output_prefix, "not exists and will be created."
		os.makedirs(os.path.dirname(args.output_prefix))
	
	
	#---------- Simulation ---------------#
	
	print "Running SimPheno..."
	covariates = None
	if args.cov != None :
		try :
			covariates = prepareCovariates(args.cov)
		except :
			print "Warning: covariates can't be prepared."
		
	## All functions below need to call from a particular class (class factory)

	# Instantiation of the Simpheno class
	sim = None
	
	try :
		sim = Simpheno(intype=args.intype)
	except :
		print sys.exc_info()
	
	try :
		sim.prepare(args.idata, args.ceff, args.output_prefix)
	except :
		print sys.exc_info()
	"""
	if args.dflag :
		# Simulate dichotomous (binary) trait
		try :
			sim.simDichotomousPhe(covariates, args.p0)
		except :
			print sys.exc_info()
	
	if args.cflag :
		# Simulate continuous (qualitative) trait
		try :
			sim.simContinuousPhe(covariates)
		except :
			print sys.exc_info()
	
	if args.sflag :
		# Simulate survival trait
		if args.weib :
			try :
				sim.simulWeib(covariates, 7e-8, 1, 0.01)
			except :
				print sys.exc_info()
		if args.gomp :
			try :
				sim.simulGomp(covariates, 7e-8, 1, 0.01)
			except :
				print sys.exc_info()
	
	
	print "Saving results..."
	try :	
		sim.saveData()
	except :
		print "Cannot save results."
		sys.exit(-1)
	"""
	print "Done!"


# Start-up code =========={{{
if __name__=="__main__":
    main()
#.........................}}}
