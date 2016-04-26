#!/usr/bin/env python

import os, sys
import rdflib


s = '''
import os, sys
sys.path.append(os.path.expanduser('~/c25ld/analysis'))
import rdflib

TrajName = '%(TRAJNAME)s'
Nbins = %(NBINS)d
Cut = %(CUT)g
Prefix = '%(PREFIX)s'

gofr = rdflib.rdf.getRDF(TrajName = TrajName, Nbins = Nbins, Cut = Cut, Prefix = Prefix,
						 frame_start = 0, frame_stop = 20000, stepfreq = 100)

'''



# Main
sJobIn = file(os.path.expanduser('~/job.template')).read()
lammpsdir = os.path.abspath(sys.argv[1])
target_dir = os.path.abspath(sys.argv[2])

fftypes = ['lj', 'wca']
cut = 20.0
nbins = 40

for fftype in fftypes:
	trajname = os.path.join(lammpsdir, 'unconstrained'+fftype, 'c25_unbiased.lammpstrj.gz')
	prefix = 'c25_%s_rdf' % fftype
	
	d1 = {'TRAJNAME' : trajname, 'NBINS' : nbins, 'CUT' : cut, 'PREFIX' : prefix}
	filename = os.path.join(target_dir, prefix + '.py')
	file(filename, 'w').write(s % d1)
	
	d2 = {'JOBNAME' : 'rdf', 'hasMAIL' : '#', 'CMD' : 'python %s' % prefix + '.py'}
	sJob = os.path.join(target_dir, prefix + '.sh') 
	file(sJob, 'w').write(sJobIn % d2)
	os.system('chmod 777 ' + sJob)


