#/usr/bin/env python

import os, sys
import numpy as np

'''
small utility script to compute scalefactor
for new box length following a NPT run crash
'''

crashlog = sys.argv[1]
lines = open(crashlog, 'r').readlines()
p = []

start = [lines.index(line) for line in lines if line.startswith('Step Atoms Press Temp')][0] + 1
stop = [lines.index(line) for line in lines if line.startswith('ERROR on proc')][0]
lines = lines[start:stop]
N = len(lines)

for i in range(N):
	l = lines[i].split()
	p.append(float(l[2]))
	
kappa_water = 46.4e-6 #compressibility of water at 298 K
delp = 1.0 - np.mean(p)
ratio = np.exp(-kappa_water * delp)
scalefactor = ratio ** (1./3.)
print scalefactor
