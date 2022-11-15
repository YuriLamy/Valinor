#!/usr/bin/python3

import pandas as pd
import numpy as np

fn = "5_49900.msat"

pf = pd.read_csv(fn,sep='\t')
nrow = len(pf["a1"])

print("МОНОГЕННЫЙ CЛУЧАЙ \n")

Ho_list = 1.0*((pf["a1"]==1) & (pf["a2"]==0) & (pf["a3"]==0)/nrow)
Ho_list += 1.0*((pf["a1"]==0) & (pf["a2"]==1) & (pf["a3"]==0)/nrow)
Ho_list += 1.0*((pf["a1"]==0) & (pf["a2"]==0) & (pf["a3"]==1)/nrow)
print ("гомозигот по гену а вcего\t",sum(Ho_list)/nrow)
print ("гетерозигот по гену а вcего\t", (nrow - sum(Ho_list))/nrow)
print("===============================================")
He_a_12 = (((1*(pf['a1'] ==1)) + (1*(pf['a2'] ==1)))==2)*1
He_a_13 = (((1*(pf['a1'] ==1)) + (1*(pf['a3'] ==1)))==2)*1
He_a_23 = (((1*(pf['a2'] ==1)) + (1*(pf['a3'] ==1)))==2)*1
print("He_a_12 = ", sum(He_a_12)/nrow)	
print("He_a_13 = ", sum(He_a_13)/nrow)	
print("He_a_23 = ", sum(He_a_23)/nrow)	

Ho_list = 1.0*((pf["b1"]==1) & (pf["b2"]==0) & (pf["b3"]==0))
Ho_list += 1.0*((pf["b1"]==0) & (pf["b2"]==1) & (pf["b3"]==0))
Ho_list += 1.0*((pf["b1"]==0) & (pf["b2"]==0) & (pf["b3"]==1))
print ("гомозигот по гену b вcего\t",sum(Ho_list))
print ("генерозигот по гену b вcего\t", nrow - sum(Ho_list))
print("===============================================")
He_b_12 = (((1*(pf['b1'] ==1)) + (1*(pf['b2'] ==1)))==2)*1
He_b_13 = (((1*(pf['b1'] ==1)) + (1*(pf['b3'] ==1)))==2)*1
He_b_23 = (((1*(pf['b2'] ==1)) + (1*(pf['b3'] ==1)))==2)*1
print("He_b_12 = ", sum(He_b_12)/nrow)	
print("He_b_13 = ", sum(He_b_13)/nrow)	
print("He_b_23 = ", sum(He_b_23)/nrow)	

Ho_list = 1.0*((pf["c1"]==1) & (pf["c2"]==0) & (pf["c3"]==0))
Ho_list += 1.0*((pf["c1"]==0) & (pf["c2"]==1) & (pf["c3"]==0))
Ho_list += 1.0*((pf["c1"]==0) & (pf["c2"]==0) & (pf["c3"]==1))
print ("гомозигот по гену c вcего\t",sum(Ho_list)/nrow)
print ("гетерозигот по гену c вcего\t", (nrow - sum(Ho_list))/nrow)
print("===============================================")
He_c_12 = (((1*(pf['c1'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
He_c_13 = (((1*(pf['c1'] ==1)) + (1*(pf['c3'] ==1)))==2)*1
He_c_23 = (((1*(pf['c2'] ==1)) + (1*(pf['c3'] ==1)))==2)*1
print("He_c_12 = ", sum(He_c_12)/nrow)	
print("He_c_13 = ", sum(He_c_13)/nrow)	
print("He_c_23 = ", sum(He_c_23)/nrow)	

print("\n ДВУХГЕННЫЙ СЛУЧАЙ\n")

a1b1 = (((1*(pf['a1'] ==1)) + (1*(pf['b1'] ==1)))==2)*1
a1b2 = (((1*(pf['a1'] ==1)) + (1*(pf['b2'] ==1)))==2)*1
a1b3 = (((1*(pf['a1'] ==1)) + (1*(pf['b3'] ==1)))==2)*1

a2b1 = (((1*(pf['a2'] ==1)) + (1*(pf['b1'] ==1)))==2)*1
a2b2 = (((1*(pf['a2'] ==1)) + (1*(pf['b2'] ==1)))==2)*1
a2b3 = (((1*(pf['a2'] ==1)) + (1*(pf['b3'] ==1)))==2)*1

a3b1 = (((1*(pf['a3'] ==1)) + (1*(pf['b1'] ==1)))==2)*1
a3b2 = (((1*(pf['a3'] ==1)) + (1*(pf['b2'] ==1)))==2)*1
a3b3 = (((1*(pf['a3'] ==1)) + (1*(pf['b3'] ==1)))==2)*1

print("a1b1 =\t", sum(a1b1)/nrow)
print("a1b2 =\t", sum(a1b2)/nrow)
print("a1b3 =\t", sum(a1b3)/nrow)
print("a2b1 =\t", sum(a2b1)/nrow)
print("a2b2 =\t", sum(a2b2)/nrow)
print("a2b3 =\t", sum(a2b3)/nrow)
print("a3b1 =\t", sum(a3b1)/nrow)
print("a3b2 =\t", sum(a3b2)/nrow)
print("a3b3 =\t", sum(a3b3)/nrow)

a1c1 = (((1*(pf['a1'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
a1c2 = (((1*(pf['a1'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
a1c3 = (((1*(pf['a1'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

a2c1 = (((1*(pf['a2'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
a2c2 = (((1*(pf['a2'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
a2c3 = (((1*(pf['a2'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

a3c1 = (((1*(pf['a3'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
a3c2 = (((1*(pf['a3'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
a3c3 = (((1*(pf['a3'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

print('\n')
print("a1c1 =\t", sum(a1c1)/nrow)
print("a1c2 =\t", sum(a1c2)/nrow)
print("a1c3 =\t", sum(a1c3)/nrow)
print("a2c1 =\t", sum(a2c1)/nrow)
print("a2c2 =\t", sum(a2c2)/nrow)
print("a2c3 =\t", sum(a2c3)/nrow)
print("a3c1 =\t", sum(a3c1)/nrow)
print("a3c2 =\t", sum(a3c2)/nrow)
print("a3c3 =\t", sum(a3c3)/nrow)

b1c1 = (((1*(pf['b1'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
b1c2 = (((1*(pf['b1'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
b1c3 = (((1*(pf['b1'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

b2c1 = (((1*(pf['b2'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
b2c2 = (((1*(pf['b2'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
b2c3 = (((1*(pf['b2'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

b3c1 = (((1*(pf['b3'] ==1)) + (1*(pf['c1'] ==1)))==2)*1
b3c2 = (((1*(pf['b3'] ==1)) + (1*(pf['c2'] ==1)))==2)*1
b3c3 = (((1*(pf['b3'] ==1)) + (1*(pf['c3'] ==1)))==2)*1

print('\n')
print("b1c1 =\t", sum(b1c1)/nrow)
print("b1c2 =\t", sum(b1c2)/nrow)
print("b1c3 =\t", sum(b1c3)/nrow)
print("b2c1 =\t", sum(b2c1)/nrow)
print("b2c2 =\t", sum(b2c2)/nrow)
print("b2c3 =\t", sum(b2c3)/nrow)
print("b3c1 =\t", sum(b3c1)/nrow)
print("b3c2 =\t", sum(b3c2)/nrow)
print("b3c3 =\t", sum(b3c3)/nrow)

a1 = sum(1*(pf['a1']==1))/nrow
a2 = sum(1*(pf['a2']==1))/nrow
a3 = sum(1*(pf['a3']==1))/nrow
b1 = sum(1*(pf['b1']==1))/nrow
b2 = sum(1*(pf['b2']==1))/nrow
b3 = sum(1*(pf['b3']==1))/nrow
c1 = sum(1*(pf['c1']==1))/nrow
c2 = sum(1*(pf['c2']==1))/nrow
c3 = sum(1*(pf['c3']==1))/nrow

print ("для каждого возможного гаплотипа и для каждой пары генов вычисляем D_ij")
DD=[
abs(sum(a1b1)/nrow-a1*b1),abs(sum(a1b2)/nrow-a1*b2),abs(sum(a1b2)/nrow-a1*b3),
abs(sum(a2b1)/nrow-a2*b1),abs(sum(a2b2)/nrow-a2*b2),abs(sum(a2b3)/nrow-a2*b3),
abs(sum(a3b1)/nrow-a3*b1),abs(sum(a3b2)/nrow-a3*b2),abs(sum(a3b2)/nrow-a3*b3),
abs(sum(a1c1)/nrow-a1*c1),abs(sum(a1b2)/nrow-a1*c2),abs(sum(a1b2)/nrow-a1*c3),
abs(sum(a2c1)/nrow-a2*c1),abs(sum(a2b2)/nrow-a2*c2),abs(sum(a2b3)/nrow-a2*c3),
abs(sum(a3c1)/nrow-a3*c1),abs(sum(a3b2)/nrow-a3*c2),abs(sum(a3b2)/nrow-a3*c3),
abs(sum(b1c1)/nrow-b1*c1),abs(sum(b1c2)/nrow-b1*c2),abs(sum(b1c2)/nrow-b1*c3),
abs(sum(b2c1)/nrow-b2*c1),abs(sum(b2c2)/nrow-b2*c2),abs(sum(b2c3)/nrow-b2*c3),
abs(sum(b3c1)/nrow-b3*c1),abs(sum(b3c2)/nrow-b3*c2),abs(sum(b3c2)/nrow-b3*c3)
]

print("среднее Dij = \t",np.average(DD))
print("максимум Dij = \t",np.max(DD))
print("минимум Dij = \t", np.min(DD))
print("варианса Dij = \t",np.var(DD))

#print("\n надо бы еще и распределение этих Д построить...")
#print(DD)
#ожидаемая частота любого аллеля любого гена, поскольку у нас модель :)
P_exp=.333
#при независимом наследовании совместная вероятность любого гаплотипа их 2 генов:
PP_exp = .333*.333