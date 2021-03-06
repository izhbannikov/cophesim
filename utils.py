#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-

from itertools import chain
import re

class Parser :

	def parse_ms(fname,diploid=0):
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
			genotypes,positions,raw=parse_ms_simwise(simlines,diploid)
			positions_all.append(positions)
			genotypes_all.append(genotypes)
			raw_all.append(raw)

		return genotypes_all,positions_all,raw_all



	def parse_ms_simwise(lines,diploid=0):
		positions=[]
		genotypes=[]

		sim_count=0

		for l in lines:
			if len(l)==0:
				continue
			if l[:3]=='pos':
				s=re.sub('positions: ','',l.rstrip())
				positions=map(float,s.split(' '))

			if l[0] in ['0','1']:
				genotypes.append(list(l.rstrip()))
		
		raw_genotypes=genotypes
			
		if diploid:
			genotypes=make_diploids(genotypes)
		

		return genotypes,positions,raw_genotypes


	def parse_genome(fname,diploid=0):
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
			genotypes,positions,raw=parse_genome_simwise(simlines,diploid)
			positions_all.append(positions)
			genotypes_all.append(genotypes)
			raw_all.append(raw)

		return genotypes_all,positions_all,raw_all
	

	def parse_genome_simwise(lines,diploid=0):

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
