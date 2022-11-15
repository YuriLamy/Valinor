#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
modelling the consequences of apomixis
each organism has mitochondrial DNA marker inherited from one of the parents 
and a set of microsats each evolving and segregating independently according simple SMM model by Kimura/Ohta 1973
"""
import sys
import random
from random import randint
import string
import time
#import matplotlib.pyplot as plt
import numpy as np

chars = ['A','C','G','T']
nmeChars = ['a','b','c','t','p','l','u','x','1','2','0','f','h','m']
'''
3 lengths of msats, SMM model 
'''
minMS = 10
maxMS = 12
msats = range(minMS,maxMS)

dirs = [-1,1]
kids = [2,3,4,5,6,7,8,9]
vegKids = 10
mitoProb = .001 	# per whole sequence!
nucProb = .0005	# per whole set of microsats!
maxNum = 10000  # no scenario at the moment, later we may do N(t) and see what happens
				# теперь Nmax размазана вокруг по Пуассону
mitoLength = 500
nucLength = 1
maxAge = 5
apoRatio = .01
lab = "apoRatio=.01"
sLength = 50000
simulationRepeats = 5
apo_immortal = 1
start_time=time.time()

###################################################################################################################

def mkName():
	res = ""
	for i in range(0,4):
		ch=random.choice(nmeChars)
		res += ch
	return res

def wConditions(root,prnt):
	fn = root+'.ctl'
	f =open(fn,'w')
	localtime = time.asctime( time.localtime(time.time()) )
	f.write("apomyx.py\nStarted at:\t%s\n"%localtime)
	f.write("max. number of organisms:\t%i\n"%maxNum)
	f.write("max.possible age:\t\t%i\n"%maxAge)
	f.write("Length of simulation:\t%i\n"%sLength)
	f.write("Probability of sexual:\t%3.2f\n"%apoRatio)
	f.write("msat mutation rate:\t%4.3f\n"%nucProb)
	f.write("mitochondrial DNA rate:\t%4.3f\n"%mitoProb)
	f.write("\n>ancestor\n")
	f.write(''.join(prnt["mito"]))
	f.write("\n")
	f.close()
	return

def observedH1(g):
	ho = 0
	he = 0
	for i in range(0,nucLength):
		if g["nuc1"][i] == g["nuc2"][i]:
			ho += 1.0
		else:
			he += 1.0
	res = he/(ho+he)
	return res
	
def observedH2(g):
	hi = 0
	hy = 0
	for i in range(0,nucLength):
		if g["nuc3"][i] == g["nuc4"][i]:
			hi += 1.0
		else:
			hy += 1.0
	res = hy/(hi+hy)
	return res
	
def observedH3(g):
	ha = 0
	hu = 0
	for i in range(0,nucLength):
		if g["nuc5"][i] == g["nuc6"][i]:
			ha += 1.0
		else:
			hu += 1.0
	res = hu/(ha+hu)
	return res

def total_observedH1(G):
	count = 0.0
	he = .0
	for g in G:
		if g["vacant"] == 0:
			count += 1.0
			he += observedH(g)
	return he/count/2

def total_observedH2(G):
	count = 0.0
	hy = .0
	for g in G:
		if g["vacant"] == 0:
			count += 1.0
			hy += observedH(g)
	return hy/count/2

def total_observedH3(G):
	count = 0.0
	hu = .0
	for g in G:
		if g["vacant"] == 0:
			count += 1.0
			hu += observedH(g)
	return hu/count/2
	
def expected_H1(G):
	haplo = list()
	loci = list()
	l = 2*len(G)
	for i in range(0,nucLength):
		tmp = list()
		haplo.append(tmp)
		for j in range(0,10):
			haplo[i].append(0.0)
	for g in G:
		for i in range(0,nucLength):
			q = g["nuc1"][i]-10
			p = g["nuc2"][i]-10
			haplo[i][q] += 1.0
			haplo[i][p] += 1.0
	for i in range(0,nucLength):
		for j in range(0,10):
			haplo[i][j] /= l
			haplo[i][j] *= haplo[i][j]
	for i in range (0,nucLength):
		loci.append(sum(haplo[i]))
	
	return 1 - sum(loci)/nucLength
	
def expected_H2(G):
	haplo = list()
	loci = list()
	l = 2*len(G)
	for i in range(0,nucLength):
		tmp = list()
		haplo.append(tmp)
		for j in range(0,10):
			haplo[i].append(0.0)
	for g in G:
		for i in range(0,nucLength):
			q = g["nuc3"][i]-10
			p = g["nuc4"][i]-10
			haplo[i][q] += 1.0
			haplo[i][p] += 1.0
	for i in range(0,nucLength):
		for j in range(0,10):
			haplo[i][j] /= l
			haplo[i][j] *= haplo[i][j]
	for i in range (0,nucLength):
		loci.append(sum(haplo[i]))
	
	return 1 - sum(loci)/nucLength
	
def expected_H3(G):
	haplo = list()
	loci = list()
	l = 2*len(G)
	for i in range(0,nucLength):
		tmp = list()
		haplo.append(tmp)
		for j in range(0,10):
			haplo[i].append(0.0)
	for g in G:
		for i in range(0,nucLength):
			q = g["nuc5"][i]-10
			p = g["nuc6"][i]-10
			haplo[i][q] += 1.0
			haplo[i][p] += 1.0
	for i in range(0,nucLength):
		for j in range(0,10):
			haplo[i][j] /= l
			haplo[i][j] *= haplo[i][j]
	for i in range (0,nucLength):
		loci.append(sum(haplo[i]))
	
	return 1 - sum(loci)/nucLength

def av_Age(G):
	sum = 0
	count = 0
	for g in G:
		if g["vacant"]==0:
			count += 1.0
			sum += g["age"]
	return sum*1.0/count

def mkname(l):
	res = ''
	for i in range(0,l):
		ch = random.choice(nmeChars)
		res += ch
	return res

def mkAncestor(l):
	res = list()
	for c in range(0,l):
		res.append(random.choice(chars))

	return res

def initGad(l1,l2):
	anc = mkAncestor(l1)
	n1 = list()
	n2 = list()
	n3 = list()
	n4 = list()
	n5 = list()
	n6 = list()
	for i in range (0,l2):
		q = random.choice(msats)
		n1.append(q)
		q = random.choice(msats)
		n2.append(q) 
		q = random.choice(msats)
		n3.append(q)
		q = random.choice(msats)
		n4.append(q) 
		q = random.choice(msats)
		n5.append(q)
		q = random.choice(msats)
		n6.append(q) 	#в случайно выбранные аллели
	sx = 0
	if random.uniform(0,1)<apoRatio:
		sx = 1
	g = {"age":0,"engaged":0,"sexy":sx,"mito":anc,"nuc1":n1,"nuc2":n2,"nuc3":n3,"nuc4":n4,"nuc5":n5,"nuc6":n6,"vacant":0}
	return g

def cloneGad(g,apo):
	res = {}
	l=len(g["mito"])
	ll = len(g["nuc1"])
	res["age"]=g["age"]
	res["engaged"]=g["engaged"]
	res["sexy"] = g["sexy"]
	res["vacant"] = 0
	res["mito"] = list()
	res["nuc1"] = list()
	res["nuc2"] = list()
	res["nuc3"] = list()
	res["nuc4"] = list()
	res["nuc5"] = list()
	res["nuc6"] = list()
	
	for i in range(0,l):
		res["mito"].append(g["mito"][i])
	for i in range(0,ll):
		res["nuc1"].append(g["nuc1"][i])
		res["nuc2"].append(g["nuc2"][i])
		res["nuc3"].append(g["nuc3"][i])
		res["nuc4"].append(g["nuc4"][i])
		res["nuc5"].append(g["nuc5"][i])
		res["nuc6"].append(g["nuc6"][i])
	if random.uniform(0,1)<apoRatio:
		res["sexy"]=1
		res["age"]=0
	return res

def mateGad(g1,g2):
	print '*',
	res = {}
	res["engaged"] = 0
	if random.uniform(0,1)<apoRatio:
		res["sexy"]=1
	res["age"]=0
	res["vacant"] = 0
	res["mito"] = list()
	res["nuc1"] = list()
	res["nuc2"] = list()
	res["nuc3"] = list()
	res["nuc4"] = list()
	res["nuc5"] = list()
	res["nuc6"] = list()
	#randomly choose the "mom" to donate mitochondrial DNA
	if random.uniform(0,1)<.5:
		res["mito"] += g1["mito"]
	else:
		res["mito"] += g2["mito"]
	#microsat markers segregate independently thus the parent chromosome/source is chosen for each one
	l = len(g1["nuc1"])
	for i in range(0,l):
		if random.uniform(0,1)<.5:
			res["nuc1"].append(g1["nuc1"][i])
		else:
			res["nuc1"].append(g1["nuc2"][i])
		if random.uniform(0,1)<.5:
			res["nuc2"].append(g2["nuc1"][i])
		else:
			res["nuc2"].append(g2["nuc2"][i])
			
		if random.uniform(0,1)<.5:
			res["nuc3"].append(g1["nuc3"][i])
		else:
			res["nuc3"].append(g1["nuc4"][i])
		if random.uniform(0,1)<.5:
			res["nuc4"].append(g2["nuc3"][i])
		else:
			res["nuc4"].append(g2["nuc4"][i])
			
		if random.uniform(0,1)<.5:
			res["nuc5"].append(g1["nuc5"][i])
		else:
			res["nuc5"].append(g1["nuc6"][i])
		if random.uniform(0,1)<.5:
			res["nuc6"].append(g2["nuc5"][i])
		else:
			res["nuc6"].append(g2["nuc6"][i])
	return res

def mitoMute(g):
	res = g.copy()
	l = len(res["mito"])
	x = random.randint(0,l-1)
	s = res["mito"]
	q = s[x]
	qq = q
	while(qq==q):
		qq=random.choice(chars)
	#print x,q,qq
	res["mito"][x]=qq
	return(res)

def msatMute(g):
	res = g.copy()
	l = len(res["nuc1"])
	x = random.randint(0,l-1)
	#x=0
	if random.uniform(0,1)<.5:
		pointer = 0
	else:
		pointer = 1
	if pointer == 0:	
		s = res["nuc1"]
	else:
		s = res["nuc2"]
	q = s[x]
	qq = q
	while qq==q:
		k = random.choice(dirs)
		qq += k
		if qq<minMS :
			qq = minMS+1
		if qq>maxMS :
			qq = maxMS-1
	#print x,q,qq
	if pointer == 0:
		res["nuc1"][x]=qq
	else:
		res["nuc2"][x]=qq
		
	res = g.copy()
	l = len(res["nuc3"])
	x = random.randint(0,l-1)
	#x=0
	if random.uniform(0,1)<.5:
		pointer = 0
	else:
		pointer = 1
	if pointer == 0:	
		s = res["nuc3"]
	else:
		s = res["nuc4"]
	q = s[x]
	qq = q
	while qq==q:
		k = random.choice(dirs)
		qq += k
		if qq<minMS :
			qq = minMS+1
		if qq>maxMS :
			qq = maxMS-1
	#print x,q,qq
	if pointer == 0:
		res["nuc3"][x]=qq
	else:
		res["nuc4"][x]=qq
		
	res = g.copy()
	l = len(res["nuc5"])
	x = random.randint(0,l-1)
	#x=0
	if random.uniform(0,1)<.5:
		pointer = 0
	else:
		pointer = 1
	if pointer == 0:	
		s = res["nuc5"]
	else:
		s = res["nuc6"]
	q = s[x]
	qq = q
	while qq==q:
		k = random.choice(dirs)
		qq += k
		if qq<minMS :
			qq = minMS+1
		if qq>maxMS :
			qq = maxMS-1
	#print x,q,qq
	if pointer == 0:
		res["nuc5"][x]=qq
	else:
		res["nuc6"][x]=qq
	return(res)


def getOlder(g,crit):
	res = g
	if res["age"] >= maxAge:
		res["vacant"] = 1
	else:
		res["age"] +=1
		if crit !=0 and g["sexy"]==1:
			res["age"] -= 1
	return(res)


def main():
	
	simulationRepeats = 3
	froot = mkName()
	for step in range(0,simulationRepeats):
		frt = froot + "_"+str(step)+'_'
		#fms = open(frt+".msat",'w')
		df = open(frt+".data","w")
		df.write("Time\tN\tHe\tHo\tHi\tHy\tHa\tHu\tavAge\tFi\n")
		gad = initGad(mitoLength,nucLength)
		wConditions(frt,gad)
		G0 = list()

		for i in range(0,1000):
			G0.append(initGad(mitoLength,nucLength))
		for g in G0:
			for j in range(0,1):
				tmp = mitoMute(g)
				g = tmp.copy()
    			for j in range(0,50):
				tmp = msatMute(g)
				g = tmp.copy()
			g["age"]=random.randint(0,maxAge)
		for tt in range(0,sLength):
			G1 = list()
			for g in G0:
				tmp = getOlder(g,apo_immortal)
				if tmp["vacant"] == 0:
					if(random.uniform(0,1)<nucProb):
						tm = msatMute(tmp)
						tmp = tm.copy()
					if(random.uniform(0,1)<mitoProb):
						tm = mitoMute(tmp)
						tmp =  tm.copy()
					G1.append(tmp)
					#если вегетативно размножаться, то тут и размножимся
					if tmp["sexy"]==0:
						for i in range(0,vegKids):
							tm = cloneGad(tmp,apoRatio)
							tm["age"]=0
							G1.append(tm)
					else:
						g["sexy"] = 0
						found = 0
						# поищем себе пару
						for q in G0:
							if q["sexy"] == 1:
								ttmp = q.copy()
								found == 1
								#q["sexy"] = 0
								q["vacant"] = 1 #из кандидатов в размножение оно уже ушло
								break
						if found == 1:
							nKids = random.choice(kids)
							for j in range(0,nKids):
								G1.append(mateGad(tmp,ttmp))
							tmp=getOlder(ttmp)
							if tmp["vacant"]==0:
								G1.append(tmp)
			nG = len(G1)*1.0
			G0 = list()
			devi = np.random.poisson(maxNum)
			if nG > devi:
				G0 = random.sample(G1,maxNum)
			else:
				G0 = G1
			
#			if tt % 50000 == 0:
#				Ho = total_observedH(G0)
#				print tt
#				He = expected_H(G0)
#				df.write( "%i\t%i\t%8.7f\t%5.4f\t%5.4f\t%5.4f\n"%(tt,len(G0),He,Ho,av_Age(G0),(He-Ho)/He))
			if tt % 100 == 0:
				fms = open(frt+str(tt)+".msat",'w')
				nms = np.random.choice(G0,1000,replace=False)
				count = 0
				for q in nms:
					fms.write("%i\t%i\t"%(step,count))
					count += 1
					
					for j in range(10,13):
						if j == q["nuc1"][0] or j == q["nuc2"][0]:
							fms.write("1\t")
							ind = 1
						else:
							fms.write("0\t")
					for j in range(10,13):
						if j == q["nuc3"][0] or j == q["nuc4"][0]:
							fms.write("1\t")
							ind = 1
						else:
							fms.write("0\t")
					for j in range(10,13):
						if j == q["nuc5"][0] or j == q["nuc6"][0]:
							fms.write("1\t")
							ind = 1
						else:
							fms.write("0\t")
					fms.write("\n")
				fms.close()
				'''
				fms = open(frt+str(tt)+".msat",'w')
				nms = np.random.choice(G0,100,replace=False)
				count = 0
				for q in nms:
					fms.write("%i\t"%count)
					count += 1
					
					for j in range(10,40):
						if j == q["nuc1"][0] or j == q["nuc2"][0]:
							fms.write("1\t")
							ind = 1
						else:
							fms.write("0\t")
					fms.write("\n")
				fms.close()
				'''
				he=0
				ho=0
				hi=0
				hy=0
				ha=0
				hu=0
				h = list()
				for g in G0:
					h.append(g["nuc1"][0])
					h.append(g["nuc2"][0])
					if(g["nuc1"][0] == g["nuc2"][0]):
						ho += 1
					else:
						he += 1
					h.append(g["nuc3"][0])
					h.append(g["nuc4"][0])
					if(g["nuc3"][0] == g["nuc4"][0]):
						hi += 1
					else:
						hy += 1
					h.append(g["nuc5"][0])
					h.append(g["nuc6"][0])
					if(g["nuc5"][0] == g["nuc6"][0]):
						ha += 1
					else:
						hu += 1
				hh = set(h)
			print("%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s"%(tt,ho,he,hi,hy,ha,hu,hh.__len__(),lab))
		df.close()
		mito = open(frt+".fas",'w')
		ll=len(G0)
		print("simulation N %i"%step)
		nms = np.random.choice(G0,1000,replace=False)
		i = 0
		for q in nms:
			mito.write(">Seq"+str(i)+"\n")
			i += 1
			sequence = ''.join(q["mito"])
			mito.write(''.join(q["mito"]))
			mito.write("\n")
		mito.close()
		#fms.close()

	print "execution time:\t%4.3f seconds"%(time.time()-start_time)
	
	return

if __name__=='__main__':
    main()



# нарисовать график на R:
#pl<-ggplot(c,aes(Time))+geom_point(aes(y=Fi),color="blue",size=.4,alpha=.5)+geom_point(aes(y=Fi1),color="magenta",size=.2)+geom_smooth(aes(y=Fi))+geom_smooth(aes(y=Fi1),color="magenta")

