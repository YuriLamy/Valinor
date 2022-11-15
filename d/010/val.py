#!/usr/bin/python3

import os,sys,re
import pandas as pd
import argparse
import numpy as np
#from Bio import Phylo
#from Bio.Phylo.TreeConstruction import DistanceCalculator
#from Bio.Phylo.TreeConstruction import DistanceMatrix
#from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def createParser ():
  parser = argparse.ArgumentParser()
  parser.add_argument ('-n', '--name', nargs='?', type = str, default='data.tab', help = 'Name of file with microsats')
  return parser

def main():
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	infile = namespace.name
	dta = pd.read_csv(infile, sep = '\t', header = 0)
	nms = pd.read_csv(nmes, sep = '\t', header = 0)
	who = dict(nms.values)
	lnumber = len(dta.s1)
	
	a = (dta.c1>0)*4(dta.c2>0)*2(dta.c3>0)
	
	s1 = sum(a & 0b100)/4
	s2 = sum(a & 0b10)/2
	s3 = sum(a & 0b1)

	print('n, s1, s2, s3')
	
if __name__=="__main__":
	main()
