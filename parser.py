#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from itertools import chain
import re
from plinkio import plinkfile

"""
class Adapter :
	# Adapter for different types of input data
	plinkParser = Parser()
	msParser = Parser()
	genomeParser = Parser()
	
	def parse() :
		# Parses input file depenting on their types
		
	
	def getData() :
		# Returns genotype matrix
		gmatrix = []
		return gmatrix
"""

class Parser() :

	def __init__(self):
		pass
	

	def parse_plink(self, fname) :
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
		
		return matrix

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


	

	