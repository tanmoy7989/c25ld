#!/usr/bin/env python

import os,sys

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

argstring = ''
for n in range(9):
    l = lines[start]
    argstring += l.split()[-1] + ' ' 
    start += 1


sJob1 = os.path.join(data_dir, 'SrelBennett.sh')
sJob2 = os.path.join(data_dir, 'dSreldrc.sh')

d1 = {'JOBNAME': 'SrelBennett', 'CMD': 'python SrelBennett.py %s' % argstring}
d2 = {'JOBNAME': 'dSreldrc' , 'CMD' : 'python dSreldrc.py %s' % argstring}

file(sJob1, 'w').write(sJobIn % d1)
file(sJob2, 'w').write(sJobIn % d2)

os.system('chmod 777 %s; chmod 777 %s' % (sJob1, sJob2))

os.system('cp %s %s' % ('dSreldrc.py', data_dir))
os.system('cp %s %s' % ('SrelBennett.py', data_dir))
