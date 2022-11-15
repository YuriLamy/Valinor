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

	a = (dta.a1>0)*4 + (dta.a2>0)*2 + (dta.a3>0)
	
	s1 = sum(a & 0b100000000)/4
	s2 = sum(a & 0b10000000)/2
	s3 = sum(a & 0b1000000)

	b = (dta.b1>0)*4 + (dta.b2>0)*2 + (dta.b3>0)
	
	s4 = sum(a & 0b100000)/4
	s5 = sum(a & 0b10000)/2
	s6 = sum(a & 0b1000)
	
	c = (dta.c1>0)*4 + (dta.c2>0)*2 + (dta.c3>0)
	
	s7 = sum(a & 0b100)/4
	s8 = sum(a & 0b10)/2
	s9 = sum(a & 0b1)

	hoab = sum(((a == 1) & (b == 1)) | ((a == 2) & (b == 2)) | ((a == 4) & (b == 4)))
	heab = sum(((a == 3) & (b == 3)) | ((a == 5) & (b == 5)) | ((a == 6) & (b == 6)))
	hoac = sum(((a == 1) & (c == 1)) | ((a == 2) & (c == 2)) | ((a == 4) & (c == 4)))
	heac = sum(((a == 3) & (c == 3)) | ((a == 5) & (c == 5)) | ((a == 6) & (c == 6)))
	hobc = sum(((b == 1) & (c == 1)) | ((b == 2) & (c == 2)) | ((b == 4) & (c == 4)))
	hebc = sum(((b == 3) & (c == 3)) | ((b == 5) & (c == 5)) | ((b == 6) & (c == 6)))

	print ( hoab,heab, hoac, hoab, hobc, hebc )
	#print(hoa, sep='\n')
	
if __name__=="__main__":
	main()