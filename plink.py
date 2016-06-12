#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

"""
A wrapper for Plink program.
Performs several tasks with plink, such as quality control, select snps, etc.
"""

import subprocess

def preparePlink(snpnull, snpdisease):
	#print "Making wgas.sim configuration file..."
	f = open("wgas.sim",'w')
	f.write("{} null 0.00 1.00 1.00 1.00\n".format(snpnull))
	f.write("{} disease 0.00 1.00 2.00 mult\n".format(snpdisease))
	f.close()


def runPlinkSimulation(path, ncases, ncontrols) :
	#plink --simulate wgas.sim --make-bed --out sim1 --simulate-ncases $N_cases --simulate-ncontrols $N_controls
	subprocess.call(["plink", "--simulate", "wgas.sim", "--make-bed", "--out", path, "--simulate-ncases", str(ncases), "--simulate-ncontrols", str(ncontrols)])

def runPlink(path, phenotype, output_prefix) :
	subprocess.call(["plink", "--noweb", "--bfile", path, "--pheno", phenotype, "--1", "--make-bed", "--out", output_prefix])
	
def plinkQC(bfile, selected_snps, individuals, outpath) :
	qc = "--mind 0.05 --hwe 0.001 --maf 0.05 --geno 0.05"
	if selected_snps != None :
		subprocess.call(["plink", "--noweb", "--bfile", bfile, "--extract", selected_snps, "--keep", individuals, "--make-bed", "--out", outpath])
		subprocess.call(["plink", "--noweb", "--bfile", outpath, "--nonfounders", "--mind", "0.05", "--hwe", "0.001", "--maf", "0.05", "--geno", "0.05", "--make-bed", "--out", outpath + "_qc"])