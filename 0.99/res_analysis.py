#!/usr/bin/python3

import os,sys,re
import pandas as pd
import argparse
import numpy as np

infile = 'data.tab'
nmes = 'names.tab'


def createParser ():
  parser = argparse.ArgumentParser()
  parser.add_argument ('-n', '--name', nargs='?', type = str, default='data.tab', help = 'Name of file with microsats')
  return parser
	
def main():
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	infile = namespace.name
	dta = pd.read_csv(infile, sep = '\t', header = 0)

	aHo100 = ((dta['a1'] == 1) & (dta['a2'] == 0) & (dta['a3'] == 0))
	print('Number of rows with first locus A (Ho) == [1,0,0]:', sum(aHo100))
	aHo010 = ((dta['a2'] == 1) & (dta['a1'] == 0) & (dta['a3'] == 0))
	print('Number of rows with first locus A (Ho) == [0,1,0]:', sum(aHo010))
	aHo001 = ((dta['a3'] == 1) & (dta['a1'] == 0) & (dta['a2'] == 0))
	print('Number of rows with first locus A (Ho) == [0,0,1]:', sum(aHo001))
	aHe110 = ((dta['a1'] == 1) & (dta['a2'] == 1) & (dta['a3'] == 0))
	print('Number of rows with first locus A (He) == [1,1,0]:', sum(aHe110))
	aHe011 = ((dta['a2'] == 1) & (dta['a1'] == 0) & (dta['a3'] == 1))
	print('Number of rows with first locus A (He) == [0,1,1]:', sum(aHe011))
	aHe101 = ((dta['a3'] == 1) & (dta['a1'] == 1) & (dta['a2'] == 0))
	print('Number of rows with first locus A (He) == [1,0,1]:', sum(aHe101))

	bHo100 = ((dta['b1'] == 1) & (dta['b2'] == 0) & (dta['b3'] == 0))
	print('Number of rows with first locus B (Ho) == [1,0,0]:', sum(bHo100))
	bHo010 = ((dta['b2'] == 1) & (dta['b1'] == 0) & (dta['b3'] == 0))
	print('Number of rows with first locus B (Ho) == [0,1,0]:', sum(bHo010))
	bHo001 = ((dta['b3'] == 1) & (dta['b1'] == 0) & (dta['b2'] == 0))
	print('Number of rows with first locus B (Ho) == [0,0,1]:', sum(bHo001))
	bHe110 = ((dta['b1'] == 1) & (dta['b2'] == 1) & (dta['b3'] == 0))
	print('Number of rows with first locus B (He) == [1,1,0]:', sum(bHe110))
	bHe011 = ((dta['b2'] == 1) & (dta['b1'] == 0) & (dta['b3'] == 1))
	print('Number of rows with first locus B (He) == [0,1,1]:', sum(bHe011))
	bHe101 = ((dta['b3'] == 1) & (dta['b1'] == 1) & (dta['b2'] == 0))
	print('Number of rows with first locus B (He) == [1,0,1]:', sum(bHe101))

	cHo100 = ((dta['c1'] == 1) & (dta['c2'] == 0) & (dta['c3'] == 0))
	print('Number of rows with first locus C (Ho) == [1,0,0]:', sum(cHo100))
	cHo010 = ((dta['c2'] == 1) & (dta['c1'] == 0) & (dta['c3'] == 0))
	print('Number of rows with first locus C (Ho) == [0,1,0]:', sum(cHo010))
	cHo001 = ((dta['c3'] == 1) & (dta['c1'] == 0) & (dta['c2'] == 0))
	print('Number of rows with first locus C (Ho) == [0,0,1]:', sum(cHo001))
	cHe110 = ((dta['c1'] == 1) & (dta['c2'] == 1) & (dta['c3'] == 0))
	print('Number of rows with first locus C (He) == [1,1,0]:', sum(cHe110))
	cHe011 = ((dta['c2'] == 1) & (dta['c1'] == 0) & (dta['c3'] == 1))
	print('Number of rows with first locus C (He) == [0,1,1]:', sum(cHe011))
	cHe101 = ((dta['c3'] == 1) & (dta['c1'] == 1) & (dta['c2'] == 0))
	print('Number of rows with first locus C (He) == [1,0,1]:', sum(cHe101))

	print('Sum AHo', sum(aHo100) + sum(aHo010) + sum(aHo001))
	print('Sum AHe', sum(aHe110) + sum(aHe011) + sum(aHe101))
	print('Sum BHo', sum(bHo100) + sum(bHo010) + sum(bHo001))
	print('Sum BHe', sum(bHe110) + sum(bHe011) + sum(bHe101))
	print('Sum CHo', sum(cHo100) + sum(cHo010) + sum(cHo001))
	print('Sum CHe', sum(cHe110) + sum(cHe011) + sum(cHe101))

	exit()

if __name__=='__main__':
    main()
