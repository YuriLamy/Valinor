#!/usr/bin/python3
# -*- coding: utf-8 -*-

file = open('f11l_0_45000.msat','r').readlines()
new_file = ['\t'.join(line.split('\t')[2:-1] if line.split('\t')[-1] == '\n' else line.split('\t')[2:]) for line in file]

f = open('01_45000.msat', 'w')
print('a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'c1', 'c2', 'c3', sep = '\t', file = f)
print(*new_file, sep = '\n', file = f)