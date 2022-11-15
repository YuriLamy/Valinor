#!/usr/bin/python3

import os,sys,re
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

infile = 'data.tab'
nmes = 'names.tab'

def main():
	dta = pd.read_csv(infile,sep='\t',header=0)
	nms = pd.read_csv(nmes,sep='\t',header=0)
	who = dict(nms.values)
	lnumber = len(dta.s1)
	#print(who)
	a = (dta.s1>0)*8+(dta.s2>0)*4+(dta.s3>0)*2+(dta.s4>0)
	s1 = sum(a & 0b1000)/8
	s2 = sum(a & 0b100)/4
	s3 = sum(a & 0b10)/2
	s4 = sum(a & 0b1)
	#unique
	us1 = sum(a == 0b1000)/s1
	us2 = sum(a == 0b100)/s2
	us3 = sum(a == 0b10)/s3
	us4 = sum(a == 0b1)/s4

	q4 = sum((a & 0b1111)==15)
	qs1 = q4/s1
	qs2 = q4/s2
	qs3 = q4/s3
	qs4 = q4/s4

	print("\\begin{table)")
	print("\\begin{flushright}\\textbf{Table header}")
	print("\\end{flushright}")
	print("\\begin{tabular}{l|p{2.5cm}p{2.5cm}p{2.5cm}}")
	print("1&2&3&4&5\\\\")
	print("\\hline")
	print("\\hline")
	print("{}&{}&{:.3f}&{:.3f}\\\\".format(who["s1"],s1,us1,qs1))
	print("{}&{}&{:.3f}&{:.3f}\\\\".format(who["s2"],s2,us2,qs2))
	print("{}&{}&{:.3f}&{:.3f}\\\\".format(who["s3"],s3,us3,qs3))
	print("{}&{}&{:.3f}&{:.3f}\\\\".format(who["s4"],s4,us4,qs4))
	print("\\end{tabular}")
	print("\\caption{text}")
	print("\\label{lab1}")
	print("\\end{table}")

	d12 = 1-sum(a == 0b1100)/(s1+s2)
	d13 = 1-sum(a == 0b1010)/(s1+s3)
	d14 = 1-sum(a == 0b1001)/(s1+s4)
	d23 = 1-sum(a == 0b110)/(s2+s3)
	d24 = 1-sum(a == 0b101)/(s2+s4)
	d34 = 1-sum(a == 0b11)/(s3+s4)
	ff = [[0,d12,d13,d14],
	[d12,0,d23,d24],
	[d13,d23,0,d34],
	[d14,d24,d34,0]]
	fft=[[0],[d12,0],[d13,d23,0],[d14,d24,d34,0]]
	
	ll = list(nms.Taxon)
	df = pd.DataFrame(ff,index = nms["Taxon"],columns=nms["Taxon"])
	ddf = DistanceMatrix(names=ll,matrix=fft)
	constructor = DistanceTreeConstructor()
	njtree = constructor.nj(distance_matrix=ddf)
	Phylo.write(njtree,"rep.tre",'newick')
	
	df.to_csv('dist.tab',sep='\t')

if __name__=="__main__":
	main()
