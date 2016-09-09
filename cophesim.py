#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

# A convenient wrapper for INTERSNP program

import subprocess
import sys
import argparse
import os
import ntpath

# Importing custom modules
from plink import *
from simpheno import *
from utils import *

# Starts the simulation
def main() :

	#------------Setting-up options for argument parser -------------#
	usage = "cophesim.py -i <path to plink files> -o <output prefix> [other options]"
	parser = argparse.ArgumentParser(usage=usage)
	
	# Input parameters
	parser.add_argument("-i", "--input", action="store", type=str, dest="idata", default=None, required=True, help="Path to input files")
	parser.add_argument("-o", "--output", action="store", type=str, dest="output_prefix", default=None, required=True, help="Output prefix")
	
	# Trait (outcome) type
	parser.add_argument("-d", "--dichotomous",action="store_true", dest="dflag", default=True, required=False, help="A flag for dichotomous phenotype, True by default.")
	parser.add_argument("-c", "--continuous",action="store_true", dest="cflag", default=False, required=False, help="A flag for continous phenotype, False by default.")
	parser.add_argument("-s", "--suvival", action="store_true", dest="sflag", default=False, required=False, help="A flag to simulate survival phenotype, False by default.")
	
	# Other parameters
	parser.add_argument("-cbetta", "--cbetta",action="store", type=float, dest="cbetta", default=0.2, required=False, help="Phenotype")
	parser.add_argument("-b0", "--baseline_mean",action="store", type=float, dest="baseline_mean", default="80", required=False, help="Baseline mean for continuous phenotype")
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
		covariates = prepareCovariates(args.cov)
		
	dich_trait = None
	cont_trait = None
	surv_trait = None
	
	## All functions below need to call from a particular class (class factory)
	
	# Get id of each individual:
	ids = getId(args.idata)	
	
	if args.dflag :
		# Simulate dichotomous (binary) trait
		dich_trait = simDichotomousPhe(args.idata, args.output_prefix, args.baseline_mean, 0.2, covariates, args.p0)
	
	if args.cflag :
		# Simulate continuous (qualitative) trait
		cont_trait = simContinuousPhe(args.idata, args.output_prefix, args.baseline_mean, 0.2, covariates)
	
	if args.sflag :
		# Simulate survival trait
		if args.weib :
			surv_trait = simulWeib(args.idata, args.output_prefix, args.baseline_mean, 0.2, covariates, 7e-8, 1, 0.01)
		if args.gomp :
			surv_trait = simulGomp(args.idata, args.output_prefix, args.baseline_mean, 0.2, covariates, 7e-8, 1, 0.01)
	
	
	print "Saving results..."	
	
	saveData(dich_trait, cont_trait, surv_trait, ids, args.output_prefix, covariates)
	
	
	print "Done!"


# Start-up code =========={{{
if __name__=="__main__":
    main()
#.........................}}}
