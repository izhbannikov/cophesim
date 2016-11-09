#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from itertools import chain
import re
import os
import imp

#----------- Importing plinkio ---------------#
try:
    plinkio_info = imp.find_module('plinkio')
    from plinkio import plinkfile
except ImportError:
	print "Warning: plinkio not installed. Install the plinkio, otherwise the program can not handle .bed, .bim and .fam files."



"""
This code was borrowed from phenosim (phenotype simulator) written by Gunter T at al.
http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-265 
and modified by Ilya Zhbannikov.
"""


class Parser() :

	def __init__(self):
		pass
	
	def parse_plink_ped(self, fname, diploid=False) :
		matrix = []
		positions = [[]]
		plink_file = open(fname, 'rU')
		lines = plink_file.readlines()
		plink_file.close()
		
		for i in range(len(lines)) :
			row = lines[i].strip('\n').split(' ')
			row = row[6:]
			matrix.append([])
			for j in range(len(row)-1):
				if diploid == False :
					if row[j] == row[j+1] :
						matrix[i].append(0)
					elif row[j] != row[j+1] :
						matrix[i].append(1)
		
		genotypes_all = [matrix]
		
		# Check for .map file
		map_file =  os.path.splitext(fname)[0] + ".map"
		if os.path.exists(map_file):
			map_file_handle = open(map_file, 'rU')
			lines = map_file_handle.readlines()
			for i in range(len(lines)) :
				row = lines[i].strip('\n').split('\t')
				[positions[0].append(row[-1])]
		else :
			[positions[0].append(0.0) for i in range(len(matrix[0]))]
		
		raw_all = [matrix]
		
		return genotypes_all,positions,raw_all	
		
		
	def parse_plink_bed(self, fname) :
		""" For parsing Plink-line files """
		matrix = []
		plink_file = plinkfile.open( fname )
		if not plink_file.one_locus_per_row( ):
			print( "This script requires that snps are rows and samples columns." )
			exit( 1 )
			
		for row in plink_file:
			for i in range(len(row)):
				matrix.append([])
			break
	
		for row in plink_file:
			for i in range(len(row)):
				matrix[i].append(row[i])
		genotypes_all = [matrix]
		positions = [[]]
		[positions[0].append(0.0) for i in range(len(matrix[0]))]
		raw_all = [matrix]
		return genotypes_all,positions,raw_all

	def parse_ms(self, fname,diploid=0):
		""" For parsing ms, msms and msHot output """
		f=open(fname)
		s=f.read()
		f.close()

		positions_all=[]
		genotypes_all=[]
		raw_all=[]

		simulations=s.split('//')[1:]

		for sim in simulations:
			simlines=sim.split('\n')
			genotypes,positions,raw=self.parse_ms_simwise(simlines,diploid)
			positions_all.append(positions)
			genotypes_all.append(genotypes)
			raw_all.append(raw)

		return genotypes_all,positions_all,raw_all



	def parse_ms_simwise(self, lines,diploid=0):
		positions=[]
		genotypes=[]

		sim_count=0

		for l in lines:
			if len(l)==0:
				continue
			if l[:3]=='pos':
				s=re.sub('positions: ','',l.rstrip())
				positions=map(float,s.split(' '))
			
			if l[0] in ['0','1', '2']:
				genotypes.append(list(l.rstrip()))
		
		raw_genotypes=genotypes
			
		if diploid:
			genotypes=make_diploids(genotypes)
		

		return genotypes,positions,raw_genotypes

	
	def parse_genome(self, fname,diploid=0):
		""" Parses GENOME output """
		f=open(fname)
		s=f.read()
		f.close()
		
		positions_all=[]
		genotypes_all=[]
		raw_all=[]

		simulations=s.split('GENOME')[1:]

		for sim in simulations:
			simlines=sim.split('\n')
			genotypes,positions,raw=self.parse_genome_simwise(simlines,diploid)
			positions_all.append(positions)
			genotypes_all.append(genotypes)
			raw_all.append(raw)

		return genotypes_all,positions_all,raw_all
		
	def parse_genome_simwise(self, lines,diploid=0):

		positions=[]
		genotypes=[]

		sim_count=0
		pos_next=0

		for i in xrange(len(lines)):
			l=lines[i]
			if len(l)==0:
				continue
			if l[:3]=='SNP':
				pos_next=1
			elif pos_next:
				s=l.rstrip()
				positions=map(float,s.split(' '))
				pos_next=0

			if l[:3]=='POP':
				s=re.sub('POP.*: ','',l.rstrip())
				genotypes.append(list(s))


		for i in xrange(len(positions)):
			positions[i]=float(positions[i])
	
		raw_genotypes=genotypes
		
		if diploid:
			genotypes=make_diploids(genotypes)


		return genotypes,positions,raw_genotypes


	

	