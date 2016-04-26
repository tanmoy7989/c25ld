#!/usr/bin/env python


import os, sys, pickle
import numpy as np

sJobIn = file(os.path.expanduser('~/job.template')).read()


def LD_distro(TrajName, LDCuts, Prefix, target_dir):

		s = '''
#!/usr/bin/env python
import os, sys, pickle
sys.path.append(os.path.expanduser('~/c25ld/analysis'))
import rdflib
		
TrajName = sys.argv[1]
Prefix = sys.argv[2]		
LDCuts = %(LDCuts)s
start = 0
stop = 10000
stepfreq = 100
LD = rdflib.rdf.getLocalDensity(TrajName = TrajName, LDAtomType = (1,1), LDUpperCuts = LDCuts,
						 Prefix = 'spce_LD', target_dir = os.getcwd(), 
						frame_start = start, frame_stop = stop, stepfreq = stepfreq)
'''
		
		# fill in the template
		d1 = {'LDCuts': str(LDCuts)}
		filename = os.path.join(target_dir ,Prefix+'.py')
		file(filename, 'w').write(s % d1)
		
		# write job file
		d2 = {'JOBNAME' : 'spce_LD', 'hasMAIL' : '#', 'CMD': 'python %s.py %s %s' % (Prefix, TrajName, Prefix)}
		sJob = os.path.join(target_dir, Prefix + '.sh')
		file(sJob, 'w').write(sJobIn % d2)
		os.system('chmod 777 %s' % sJob)
		
			

### Main			
									
target_dir = os.path.abspath('../data/analysis/water_spce')
TrajName = os.path.abspath('../data/lammpstraj/water_spce/spce.lammpstrj.gz')

# Solvation shells in SPC/E water are: {3.0, 5.0, 7.0, 8.0}
# Also minimum distance of approach = 2*1.4 = 2.8
LDCuts = [3., 5., 7., 8., 10., 12., 15., 18., 20., 22., 25., 30.]
Prefix = 'spce_LD'			
LD_distro(TrajName, LDCuts, Prefix, target_dir)		
