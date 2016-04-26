#!/usr/bin/env python

import os, sys, pickle
import numpy as np

sJobIn = file(os.path.expanduser('~/job.template')).read()

		
def FirstShellWaters(TrajName, LDCuts, FirstShellCut, Prefix, target_dir):

	s = '''
#!/usr/bin/env python
import os, sys, pickle
sys.path.append(os.path.expanduser('~/c25ld/analysis'))
import rdflib
	
TrajName = sys.argv[1]
Prefix = sys.argv[2]	
start = 0
stop = 20000
stepfreq = 100
(Nw, LD) = rdflib.rdf.getFirstShellWater(TrajName = TrajName, LDUpperCuts = %(LDCuts)s, Prefix = Prefix, FirstShellCut = %(FirstShellCut)g,
									  frame_start = start, frame_stop = stop, stepfreq = stepfreq)
'''
	
	# fill in the template
	d1 = {'LDCuts': str(LDCuts), 'FirstShellCut': FirstShellCut}
	filename = os.path.join(target_dir, Prefix + '.py')
	file(filename, 'w').write(s%d1)
	
	# write job file
	d2 = {'JOBNAME': 'firstshell' , 'hasMAIL': '#', 'CMD': 'python %s.py %s %s' % (Prefix, TrajName, Prefix)}
	sJob = os.path.join(target_dir, Prefix + '.sh')
	file(sJob, 'w').write(sJobIn % d2)
	os.system('chmod 777 %s' % sJob)
	
		
		

### Main

fftypes = ['lj', 'wca']
lammpsdir = os.path.abspath(sys.argv[1])
target_dir = os.path.abspath(sys.argv[2])

for fftype in fftypes:
	Prefix = 'firstshell_%s' % fftype
	
	TrajName  = os.path.join(lammpsdir, 'unconstrained_%s' % fftype, 'c25_unbiased.lammpstrj.gz')
	
	## LDCuts and FirstShellCut are based on monomer-water rdf solvation shells
	## Data and plots for these rdfs are kept in ~/c25ld/data/analysis/feb15_runs/rdf
	FirstShellCut = 3.5;
	LDCuts = [5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 15.0]

	FirstShellWaters(TrajName, LDCuts, FirstShellCut, Prefix, target_dir)
