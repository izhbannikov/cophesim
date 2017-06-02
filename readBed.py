#!/usr/bin/env vpython
#   -*- coding: utf-8 -*-


import cPickle

def readBed(filepath):
    with open(filepath, 'r') as f:
        f.seek(0)
        data = cPickle.load(f)
        return data

def writeBed(models, filepath):
    with open(filepath, 'w') as f:
        cPickle.dump(models, f)
        
        
data = readBed("/Users/ilya/Projects/cophesim/sim.plink.bed")