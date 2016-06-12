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
	
	# Setting-up options for option parser:
	usage = "cophesim.py -in <path to plink files> -out <output prefix> [other options]"
	parser = argparse.ArgumentParser(usage=usage)
	parser.add_argument("-in", "--input", action="store", type=str, dest="plink_path", default=None, required=True, help="Path to plink files")
	parser.add_argument("-out", "--output", action="store", type=str, dest="output_prefix", default=None, required=True, help="Output prefix")
	parser.add_argument("-ncase", "--ncase",action="store", type=int, dest="ncase", default=1000, required=False, help="Number of cases") 
	parser.add_argument("-ncont", "--ncont",action="store", type=int, dest="ncontrol", default=1000, required=False, help="Number of controls") 
	parser.add_argument("-snpn", "--snpn",action="store", dest="snp_null", default=100, required=False, help="Number of null SNPs") 
	parser.add_argument("-snpd", "--snpd",action="store", type=int, dest="snp_disease", default="10", required=False, help="Number of disease SNPs")
	parser.add_argument("-cbetta", "--cbetta",action="store", type=float, dest="cbetta", default=0.2, required=False, help="Phenotype")
	parser.add_argument("-b0", "--baseline_mean",action="store", type=float, dest="baseline_mean", default="80", required=False, help="Baseline mean for continuous phenotype")
	args = parser.parse_args()
	#-------------------------#
	
	
	
	# Simulation genotypes with Plink.
	# For help visit: http://pngu.mgh.harvard.edu/~purcell/plink/simulate.shtml
	preparePlink(args.snp_null, args.snp_disease)
	runPlinkSimulation(args.plink_path, args.ncase, args.ncontrol)
	
	print "Running SimPheno"
	cont_trait = simContinuousPhe(args.plink_path, args.plink_path + ".simfreq", args.baseline_mean, args.cbetta)
	dich_trait = simDichotomousPhe(args.plink_path)
	ids = getId(args.plink_path)
	print "Saving results..."
	saveData(cont_trait, dich_trait, ids, args.output_prefix)
	print "Done!"


# Start-up code =========={{{
if __name__=="__main__":
    main()
#.........................}}}
