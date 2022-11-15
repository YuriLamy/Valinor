#!/usr/bin/python3

import os,sys,re
import pandas as pd
import argparse
import numpy as np

infile = 'data.tab'
nmes = 'names.tab'


def createParser ():
  parser = argparse.ArgumentParser()
  parser.add_argument ('-n', '--name', nargs='?', default='data.tab')
  return parser
	
def main():
	parser = createParser()
	namespace = parser.parse_args(sys.argv[1:])
	infile = namespace.name
	dta = pd.read_csv(infile, sep = '\t', header = 0)
	bb = ((dta['a1'] == 1) & (dta['a2'] == 0) & (dta['a3'] == 0))
	#print(bb)
	print('Number of rows with first locus == [1,0,0]:', sum(bb))
	exit()

	"""
	nms = pd.read_csv(nmes,sep='\t',header=0)
	lnumber = len(dta["s1"])
	dd = pd.DataFrame(0,index = np.arange(lnumber),columns=nms["Name"])
	ddist = pd.DataFrame(0,index = nms["Name"],columns=nms["Name"])
	columns = list(dd)
	for z in columns:
		dd[z] = 1*(dta[z]>0)
		
	a = 1*((dta['s1']>0) & (dta['s2']==0))
	for z in nms["Name"]:
		for u in nms["Name"]:
			ddist[z][u] =1- 2*sum(1*((dta[z]>0) & (dta[u]>0)))/(sum(1*(dta[z]>0))+sum(1*(dta[u]>0)))
	print(ddist)
	bb = 1*(((dta['s1']>0) & (dta['s2']==0)) & (dta['s3']>0) & (dta['s4']>0))
	print('только s2==0\t', sum(bb))
	bb = 1*(((dta['s2']>0) & (dta['s1']==0)) & (dta['s3']>0) & (dta['s4']>0))
	print('только s1==0\t'+str(sum(bb)))
	bb = 1*(((dta['s1']>0) & (dta['s3']==0)) & (dta['s2']>0) & (dta['s4']>0))
	print('только s3==0\t'+str(sum(bb)))
	bb = 1*(((dta['s1']>0) & (dta['s4']==0)) & (dta['s2']>0) & (dta['s3']>0))
	print('только s4==0\t'+str(sum(bb)))
	print()
	for z in columns:
		print(z+'\t'+str(sum(dd[z])))		
	
	summas = [0,0,0,0]	
	for i in range(0,lnumber):
		s=0
		for z in columns:
			s+=dd[z][i]
		summas[s-1] += 1
	print(summas)		"""
if __name__=='__main__':
    main()
