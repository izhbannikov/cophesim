#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

import subprocess
import argparse
import os
import ntpath
import sys
import traceback

# Importing custom modules
from simpheno import *
from utils import *

version = "1.4.1"

# Starts the simulation
def main() :
    nsymb = 60
    
    #------------ Checking for dependencies --------------#
    try:
        plinkio_info = imp.find_module('plinkio')
    except ImportError:
        print "Warning: plinkio not installed. Install the plinkio, otherwise the program can not handle .bed, .bim and .fam files."


    #------------Setting-up options for argument parser -------------#
    usage = "cophesim.py -i <path to plink files> -o <output prefix> [other options]"
    parser = argparse.ArgumentParser(usage=usage)
	
    parser.add_argument("-ver", "--version", action="store_true", dest="ver", default=False, required=False, help="Version number.")
    
    # Input parameters #
    parser.add_argument("-i", "--input", action="store", type=str, dest="idata", default=None, required=False, help="Path with prefix to your input file(s). OR to your original database if mode='sampling'")
    parser.add_argument("-o", "--output", action="store", type=str, dest="output_prefix", default=None, required=False, help="Output prefix.")
    parser.add_argument("-itype", action="store", type=str, dest="itype", default="plink", required=False, help="Input format: plink (for Plink, default), ms (for ms, msms, msHot), genome (for Genome). Default: 'plink'")
	
	
    # Trait (outcome, phenotype) type #
    parser.add_argument("-d", "--dichotomous",action="store_true", dest="dflag", default=True, required=False, help="A flag for dichotomous phenotype, True by default.")
    parser.add_argument("-c", "--continuous",action="store_true", dest="cflag", default=False, required=False, help="A flag for continous phenotype, False by default.")
    parser.add_argument("-s", "--suvival", action="store_true", dest="sflag", default=False, required=False, help="A flag to simulate survival phenotype, False by default.")
	
    # Effects from causal SNPs
    parser.add_argument("-ce", action="store", dest="ceff", default=None, required=False, help="A path to the file with effect of each causal SNP. Must be in format: snp_index:effect. One snp per line")
	
    # Parameters for the output format #
    parser.add_argument("-otype", action="store", dest="otype", default="plink", required=False, help="Indicates output format, default=plink. Other possible output format: blossoc (for BLOSSOC), qtdt (for QTDT), tassel (for Tassel), emmax (for EMMAX).")
	
    # Miscellaneous parameters #
    #parser.add_argument("-hh", action="store", type=float, dest="h", default=0.8, required=False, help="TODO")
    parser.add_argument("-alpha", action="store", type=float, dest="alpha", default=0.2138, required=False, help="An 'alpha' parameter for inverse probability equation for the Gompertz hazard (see Bender at al., Generating survival times to simulate Cox proportional hazards models), 2005.")
    parser.add_argument("-epi", action="store", type=str, dest="epifile", default=None, required=False, help="File with interacting SNPs. One pair per line. Format: snp1_index,snp2_index,effect")
    parser.add_argument("-weib", action="store_true", dest="weib", default=True, required=False, help="A flag to use Weibull distribution for survival phenotype. True by default.")
    parser.add_argument("-gomp", action="store_true", dest="gomp", default=False, required=False, help="A flag to use Gompertz distribution for survival phenotype. False by default.")
    
    # To simulate LD #
    parser.add_argument("-LD", action="store", type=str, dest="ldfile", default=None, required=False, help="File with collinear SNPs. One pair per line. Format: snp1_index,snp2_index,correlation_coeff([-1,1])")
    
    # Mode #
    parser.add_argument("-mode", action="store", type=str, dest="smode", default="default", required=False, help="Sets the mode. Can be 'default', 'asis' or 'sampling'.")
	
    #parser.add_argument("-nind", action="store", type=str, dest="smode", default="default", required=False, help="Sets the mode. Can be 'default', 'asis' or 'sampling'.")
	
	
    #--------- Checking the input arguments ----------#
    args = parser.parse_args()
    
    if args.ver :
        print "Version:", version
        return
	
    if not os.path.exists(os.path.dirname(args.idata)) and args.itype != "plink":
        print "Directory", args.indata, "not exists. Aborting."
        sys.exit(-1)
	
    if not os.path.exists(os.path.dirname(args.output_prefix)):
        print "Directory", args.output_prefix, "not exists and will be created."
        if len(os.path.dirname(args.output_prefix)) != 0 :
            os.makedirs(os.path.dirname(args.output_prefix))
	
	
    #---------- Simulation ---------------#
	
    print "Running cophesim..."
		
    ## All functions below need to call from a particular class (class factory)

    # Instantiation of the Simpheno class
    sim = None
    sim = Simpheno(inargs=args)
    
    if args.smode == "sampling" :
        try :
            print "Make artificial genetic data by sampling real dataset..."
            sim.prepare(args.idata, args.ceff, args.output_prefix)
            sim.simulateBySampling()
            print "Done"
        except :
		    print sys.exc_info()
    elif args.smode == "default" :
        # Simulate phenotype
        try :
            
            print "Phenotype simulation..."
            sim.prepare(args.idata, args.ceff, args.output_prefix)
            
            if args.dflag :
                sim.simDichotomousPhe()
            if args.cflag :
                sim.simContinuousPhe()
            if args.sflag :
                if args.weib :
                    sim.simulWeib(7e-8, 1, 0.01)
                if args.gomp :
                    sim.simulGomp(7e-8, 1, 0.01)
        except :
            print "Exception in user code:"
            print '-'*nsymb
            traceback.print_exc(file=sys.stdout)
            print '-'*nsymb
			
    #-----------------Saving results and summary statistics--------------------#              
    print "Saving results..."
    
    try :	
        sim.saveData()
    except :
        print "Exception in user code:"
        print '-'*nsymb
        traceback.print_exc(file=sys.stdout)
        print '-'*nsymb
		
    print "Summary statistics..."
    try :
        sim.summary()
    except :
        print "Exception in user code:"
        print '-'*nsymb
        traceback.print_exc(file=sys.stdout)
        print '-'*nsymb

    print "Done!"


# Start-up code =========={{{
if __name__=="__main__":
    main()
#.........................}}}
