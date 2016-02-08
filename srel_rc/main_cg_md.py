#!/usr/bin/env python

import os, sys
import numpy as np

LDCuts = np.arange(4.5, 10.0, 0.1)
ArgFile = '/home/cask0/home/tsanyal/srel_rc/cmdlinargs.txt'
token = '>>>'
sJobIn = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s
date
%(CMD)s
'''

data_dir = os.path.abspath(sys.argv[1])
lines = open(ArgFile, 'r').readlines()
start = [lines.index(line) for line in lines if line.startswith(token) and line.split()[-1] == data_dir][0]
start += 1

cgArgString0 = ''
mdArgString0 = ''

ind = start
for n in range(8):
    if n<7:    cgArgString0 += lines[ind].split()[-1] + ' '
    mdArgString0 += lines[ind].split()[-1] + ' ' 
    ind += 1

for n, LDCut in enumerate(LDCuts):
    sJob = os.path.join(data_dir, 'cg_md_%d.sh' % n)
    s1 = '%s %g %d' % (cgArgString0, LDCut, n)
    s2 = '%s %g %d' % (mdArgString0, LDCut, n)
    cmdstring = 'python cg.py %s\npython md.py %s\n' % (s1,s2)
    d = {'JOBNAME': 'cg_md', 'CMD': cmdstring}  
    file(sJob, 'w').write(sJobIn % d)
    os.system('chmod 777 %s' % sJob)

os.system('cp %s %s' % ('cg.py', data_dir))
os.system('cp %s %s' % ('md.py', data_dir))
    
