#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

import imp
from math import log10, sqrt, log, exp
import os
from random import gauss, uniform, expovariate

#----------- Importing plinkio ---------------#
try:
    plinkio_info = imp.find_module('plinkio')
    from plinkio import plinkfile
except ImportError:
	print "Warning: plinkio not installed. Install the plinkio, otherwise the program can not handle .bed, .bim and .fam files."


#----------- Importing Pandas ---------------#
try:
    pandas_info = imp.find_module('numpy')
    import numpy
except ImportError:
    print "Warning: pandas not installed. Install the pandas."


#----------- Importing custom modules -----------------#
from exceptions import *
from parser import Parser
from writer import Writer


class Simpheno():
    def __init__(self, inargs):
        #---------------------------------------------------------#
        self.description = "Main class for phenotype simulation"
        self.author = "Ilya Y. Zhbannikov"
        #---------------------------------------------------------#
        # Initialising class members #
        self.inType = inargs.itype
        self.oType = inargs.otype
        # Lists to save the results
        self.dtrait = None
        self.ctrait = None
        self.strait = None
        self.indiv_id = None
        # Paths #
        self.datapath = None
        self.cefile = None
        self.output_prefix = None
        # Genetics #
        self.genoMatrix = None
        self.alleleFreq = None
        self.snpeff = dict()
        self.epieff = dict()
        self.ldeff = dict()
        # Misc #
        #self.h = inargs.h
        self.alpha = inargs.alpha
        self.epifile = inargs.epifile
        self.ldfile = inargs.ldfile
        
	
    #----------------- Methods ------------------#
    def mean(self, numbers):
        """Calculates average value of provided numbers"""
        return float(sum(numbers)) / max(len(numbers), 1)

	
    def prepareCE(self):
        """ Prepares a list of effects from causal SNPs """
		
        eff = dict()
		
        # Read cefile
        f = open(self.cefile)
        lines = f.readlines()
        f.close()
		
        for l in lines :
            splits = l.split(':')
            eff[int(splits[0])] = float(splits[1].replace('\n',''))
		
        self.snpeff = eff

    def prepareLD(self):
        """ Prepares a list of ld between SNPs """
        ld = dict()
		
        # Read cefile
        f = open(self.ldfile)
        lines = f.readlines()
        f.close()
		
        for l in lines :
            splits = l.split(',')
            if len(splits) == 3 :
                snp1 = int(splits[0].strip())
                snp2 = int(splits[1].strip())
                value = float(splits[2].replace('\n','').strip())
                if snp1 not in ld :
                    ld[snp1] = dict()
                    ld[snp1][snp2] = value
                else :
                    ld[snp1][snp2] = value
		
        self.ldeff = ld
	
    def prepareEpi(self):
        """ Prepares a list of effects from causal SNPs """
        epi = dict()
        # Read cefile
        f = open(self.epifile)
        lines = f.readlines()
        f.close()
		
        for l in lines :
            splits = l.split(',')
            if len(splits) == 3 :
                snp1 = int(splits[0].strip())
                snp2 = int(splits[1].strip())
                effect = float(splits[2].replace('\n','').strip())
                if snp1 not in epi :
                    epi[snp1] = dict()
                    epi[snp1][snp2] = effect
                else :
                    epi[snp1][snp2] = effect
		
        self.epieff = epi
		
    def prepareMatrix( self, fname ):
        """Prepares genotype matrix"""
        p = Parser()
		
        if self.inType == "plink" :
            ext = os.path.splitext(fname)[1]
            if ext.lower() == ".ped" :
                self.genoMatrix = p.parse_plink_ped(fname) #p.parse_plink(fname)
            else :
                self.genoMatrix = p.parse_plink_bed(fname)
        elif self.inType == "ms" :
            self.genoMatrix = p.parse_ms(fname)
        elif self.inType == "genome" :
            self.genoMatrix = p.parse_genome(fname)
            
        
        """ Prepare genotype frequencies """
        freq = []
		
        for j in range(len(self.genoMatrix[0][0][0])) :
            z = 0.0
            o = 0.0
            t = 0.0
            for i in range(len(self.genoMatrix[0][0])) :
                snp = float(self.genoMatrix[0][0][i][j])
                if snp == 0 :
                    z += 1.0
                elif snp == 1 :
                    o += 1.0
                elif snp == 2 :
                    t += 1.0
            
            nrow = len(self.genoMatrix[0][0])
            fA = z/nrow + 0.5*o/nrow
            fB = t/nrow + 0.5*o/nrow
            
            if fA == 0.0 :
                fA = 1e-12
                fB = 1 - fA
            if fB == 0.0 :
                fB = 1e-12
                fA = 1 - fB
            
            freq.append(min(fA,fB))
        
        self.alleleFreq = freq	

	def getId(self):
		"""Extracts ids from provided plink-formatted file"""
		plink_file = plinkfile.open( self.datapath )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			raise plinkioReadError
		
		sample_list = plink_file.get_samples( )
		if not sample_list :
			raise plinkioGetSamplesError
	
		ids = []
		for sample in sample_list:
			#print sample.fid, sample.iid, sample.phenotype
			ids.append([sample.fid, sample.iid])
	
		self.indiv_id = ids

	def makeLD(self) :
		""" Simulate LD (multicollinearity) """
		for k in self.ldeff.keys() :
			for kk in self.ldeff[k].keys() :
				#print k,kk,self.ldeff[k][kk]
				#print len(self.genoMatrix[0][0][0])
				for i in xrange(len(self.genoMatrix[0][0])):
					#for j in range(len(self.genoMatrix[0][0][i])) :
					self.genoMatrix[0][0][i][k] = self.genoMatrix[0][0][i][kk]
					#print i
					#print i/float(len(self.genoMatrix[0][0]))
					#print self.ldeff[k][kk]
					if i/float(len(self.genoMatrix[0][0])) >= self.ldeff[k][kk] :
						break

    def simDichotomousPhe(self, cov=None) :
        """Simulates dichotomous phenotype"""
        bin_trait = []
		
        if cov != None:
            pass # TODO
        else :
            if len(self.ldeff) != 0 :
                self.makeLD()
			
            for row in self.genoMatrix[0][0] :
                sum_wij_ui = 0.0
                for g,freq,j in zip(row, self.alleleFreq, range(len(row))) : #, locus_list) :
                    g = float(g)
                    wij = 0.0
                    val = 0
                    if len(self.snpeff) == 0 :
                        wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
                        if j in self.ldeff :
                            for jj in self.ldeff[j].keys() :
                                val = g*row[jj]*self.ldeff[j][jj]*wij
						
                        sum_wij_ui += val
                    else :
                        if j in self.snpeff :
                            #print self.ldeff.keys()
                            #if j in self.ldeff :
                            #	#print self.ldeff[j].keys()
                            #	for jj in self.ldeff[j].keys() :
                            #		val += g*row[jj]*self.ldeff[j][jj]*self.snpeff[j]
                            #else :
                            #	val = g*self.snpeff[j]
                            val = g*self.snpeff[j]
                            sum_wij_ui += val
                    
                    # Interactions
                    if j in self.epieff :
                        for jj in self.epieff[j].keys() :
                            sum_wij_ui += g*row[jj]*self.epieff[j][jj]
			
                z = sum_wij_ui + gauss(0, 1)
                pr = 1/(1+exp(-z))
                bt = 1 if pr > uniform(0,1) else 0
				
                #bin_trait.append([sample.fid, sample.iid, bt])
                bin_trait.append(bt)
        self.dtrait = bin_trait

    def simContinuousPhe(self, cov=None):
        """Used to simulate continuous phenotype"""
        cont_trait = []
        
        if cov != None:
            pass # TODO
        else :
            if len(self.ldeff) != 0 :
                self.makeLD()
                
            for row in self.genoMatrix[0][0] :
                sum_wij_ui = 0.0
                for g,freq,j in zip(row, self.alleleFreq, range(len(row))) :
                    g = float(g)
                    wij = 0.0
                    val = 0
                    if len(self.snpeff) == 0 :
                        wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
						#if j in self.ldeff :
						#	for jj in self.ldeff[j].keys() :
						#		val = g*row[jj]*self.ldeff[j][jj]*wij
						
                        sum_wij_ui += val
                    else :
                        if j in self.snpeff :
                            #print self.ldeff.keys()
                            #if j in self.ldeff :
                            #	#print self.ldeff[j].keys()
                            #	for jj in self.ldeff[j].keys() :
                            #		val += g*row[jj]*self.ldeff[j][jj]*self.snpeff[j]
                            #else :
                            #	val = g*self.snpeff[j]
                            val = g*self.snpeff[j]
                            sum_wij_ui += val
							
                    # Interactions
                    if j in self.epieff :
                        for jj in self.epieff[j].keys() :
                            sum_wij_ui += g*row[jj]*self.epieff[j][jj]
                
                z = sum_wij_ui + gauss(0, 1)	
				
                cont_trait.append(z)
			
        self.ctrait = cont_trait

	
	def prepareCovariates(cov):
		covar = []
		mean_cov = prepareMeanCov(cov)
		for i in range(n):
			row = []
			for m in mean_cov :
				row.append(gauss(m, sqrt(m/10)))
			covar.append(row)
		return covar
	
	
    def simulWeib(self, lambd, rho, rateC, cov=None):
        """Simulates survival time with Weibul distribution
           lamb = scale parameter in h0()
           rho = shape parameter in h0()
           rateC = rate parameter of the exponential distribution of C
        """
        surv_trait = []
		
        if len(self.ldeff) != 0 :
            self.makeLD()
		
        for row in self.genoMatrix[0][0]:
            sum_wij_ui = 0.0
            for g,freq,j in zip(row,self.alleleFreq, range(len(row))) :
                g = float(g)
                wij = 0.0
                if len(self.snpeff) == 0 :
                    wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
                    sum_wij_ui += g*wij
                else :
                    if j in self.snpeff :
                        sum_wij_ui += g*self.snpeff[j]
						
                # Interactions
                if j in self.epieff :
                    for jj in self.epieff[j].keys() :
                        sum_wij_ui += g*row[jj]*self.epieff[j][jj]
				
            v = uniform(0,1)
			
            Tlat = pow((-1*log(v) / (lambd * exp(sum_wij_ui))),(1/rho))
            # censoring times
            C = expovariate(rateC)
            # follow-up times and event indicators
            time = min(Tlat, C)
            status = 0 if Tlat >= C else 1
            surv_trait.append([time, status])
            #surv_trait.append(time)
        self.strait = surv_trait
  	
  	
    def simulGomp(self, lambd, rho, rateC, cov=None):
        """Simulates Gompertz latent event times 
           lamb = scale parameter in h0()
           rho = shape parameter in h0()
           rateC = rate parameter of the exponential distribution of C
        """
        surv_trait = []
		
        if len(self.ldeff) != 0 :
            self.makeLD()
        
        for row in self.genoMatrix[0][0]:
            sum_wij_ui = 0.0
            for g,freq,j in zip(row,self.alleleFreq, range(len(row))) :
                g = float(g)
                wij = 0.0
                if len(self.snpeff) == 0 :
                    wij = (g - 2.0*freq) / sqrt(2.0*freq*(1.0 - freq))
                    sum_wij_ui += g*wij
                else :
                    if j in self.snpeff :
                        sum_wij_ui += g*self.snpeff[j]
				
                # Interactions
                if j in self.epieff :
                    for jj in self.epieff[j].keys() :
                        sum_wij_ui += g*row[jj]*self.epieff[j][jj]
				
            v = uniform(0,1)
            Tlat = 1/self.alpha * log(1 - (self.alpha*log(v)) / (lambd * exp(sum_wij_ui)) )
            # censoring times
            C = expovariate(rateC)
            # follow-up times and event indicators
            time = min(Tlat, C)
            status = 0 if Tlat >= C else 1
            
            surv_trait.append([time, status]) # surv_trait.append(time)
        
        self.strait = surv_trait
  		
    def prepare(self, path, cefile, outprefix) :
        self.datapath = path
        self.cefile = cefile
        self.output_prefix = outprefix
		
        # Prepare SNP effects from file
        if self.cefile != None :
            self.prepareCE()
			
        # Prepare interaction effects from provided file
        if self.epifile != None :
            self.prepareEpi()
		
        # Prepare linkage disequlibriums 
        if self.ldfile != None :
            self.prepareLD()
	
        self.prepareMatrix(path)
		
    def saveData(self):
        if self.dtrait != None :
            # Saving dichotomous phenotype #
            ofname = self.output_prefix + "_pheno_bin.txt"
            f = open(ofname, "w")
            f.write("id sex pheno\n")
            iid = []
            [iid.append(i) for i in range(len(self.dtrait))]
            for id, d in zip(iid, self.dtrait) :
                f.write( "%s 0 %s\n" % (str(id), str(d)) ) 
            f.close()
            # Formatted writing #
            self.writeInSpecificFormat(ofname, self.dtrait)
	
        if self.ctrait != None :
            # Saving continuous phenotype #
            ofname = self.output_prefix + "_pheno_cont.txt"
            f = open(ofname, "w")
            f.write("id sex pheno\n")
            iid = []
            [iid.append(i) for i in range(len(self.ctrait))]
            for id, c in zip(iid, self.ctrait) :
                f.write( "%s 0 %s\n" % (str(id), str(c)) ) 
            f.close()

            # Formatted writing #
            self.writeInSpecificFormat(ofname, self.ctrait)

        if self.strait != None :
            # Saving continuous phenotype #
            ofname = self.output_prefix + "_pheno_surv.txt"
            f = open(ofname, "w")
            f.write("id sex age case\n")
            iid = []
            [iid.append(i) for i in range(len(self.strait))]
            for id, s in zip(iid, self.strait) :
                f.write( "%s 0 %s %s\n" % (str(id), str(s[0]), str(s[1])) ) 
            f.close()
			
            # Formatted writing #
            trait = []
            [trait.append(self.strait[i][0]) for i in range(len(self.strait))]
            self.writeInSpecificFormat(ofname, trait)


        if self.indiv_id != None :
            # Saving fid and id #
            ofname = self.output_prefix + "_id.txt"
            f = open(ofname, "w")
            for i in self.indiv_id :
                f.write('\t'.join(map(str,i)))
                f.write('\n')
            f.close()
			
	
		
		
    def writeInSpecificFormat(self, ofname, trait) :
        w = Writer()
        if self.oType == "emmax" :
            w.convert2emma(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
        elif self.oType == "plink" :
            w.convert2plink(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname, diploid=1)
        elif self.oType == "blossoc" :
            w.convert2blossoc(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
        elif self.oType == "qtdt" :
            w.convert2qtdt(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
        elif self.oType == "tassel" :
            w.convert2tassel(self.genoMatrix[0][0],self.genoMatrix[1][0],trait,ofname)
			
    def summary(self) :
        """ Generates summary statistics """
        ofname = self.output_prefix + "_summary.txt"
        f = open(ofname, "w")
        # Provided input parameters #
        f.write("Input type: " + self.inType+'\n')
        f.write("Output type: " + self.oType+'\n')
        f.write("Input data: " + self.datapath+'\n')
        f.write("Cefile: " + str(self.cefile)+'\n')
        f.write("Epifile: " + str(self.epifile)+'\n')
        f.write("Output prefix: " + self.output_prefix+'\n')
        #f.write("h: " + str(self.h)+'\n')
        f.write("alpha: " + str(self.alpha)+'\n')
        # Results summary #
        f.write( "Number of SNPs: %s\n" % str(len(self.genoMatrix[0][0][0])) )
        f.write( "Number of individuals: %s\n" % str(len(self.genoMatrix[0][0])) )
        # For number of cases and controls #
        ncase = sum(self.dtrait)
        ncont = len(self.genoMatrix[0][0]) - ncase
        f.write( "Number of cases: %s\n" % str(ncase) )
        f.write( "Number of controls: %s\n" % str(ncont) )
        # Average for continous phenotype #
        if self.ctrait != None :
            f.write( "Average for continous phenotype: %s\n" % str(self.mean(self.ctrait)) )
        # Average for survival phenotype #
        #if self.strait != None :
        #	f.write( "Average for survival phenotype: %s" % str(self.mean(self.strait)) )
        f.close()

    def simulateBySampling(self) :
        # Step 1: genetic simulation #
        ## Data sampling ##
        nrow = len(self.genoMatrix[0][0])
        ncol = len(self.genoMatrix[0][0][0])
        new_matrix = []
        genoMatrix = [map(list, zip(*self.genoMatrix[0][0]))]
        
        for j in range(ncol) :
            new_col = numpy.random.choice(genoMatrix[0][j], size=nrow, replace=True, p=None).tolist()
            #for i in range(nrow) :
                #new_snp = numpy.random.choice(self.genoMatrix[0][0][j], size=1, replace=True, p=None)
            #new_col.append(new_snp)
            #print new_col
            new_matrix.append(new_col) 
        
        new_matrix = map(list, zip(*new_matrix))
        genotypes_all = [new_matrix]
        positions = [[]]
        [positions[0].append(0.0) for i in range(len(new_matrix[0]))]
        raw_all = [new_matrix]
        self.genoMatrix = [genotypes_all, positions, raw_all]
        
	    # Step 2: now we are ready for phenotype simulation #
	