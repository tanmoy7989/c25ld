#/usr/bin/env python

import os,sys
import numpy as np

src_dir = sys.argv[1]
target_dir = sys.argv[2]

sJobIn = file(os.path.expanduser('~/job.template')).read()
d = {}

fftypes = ['lj', 'wca']

for fftype in fftypes:
	lammpstraj = os.path.abspath(os.path.join(src_dir, 'unconstrained_'+fftype, 'c25_unbiased.lammpstrj.gz'))
	d['JOBNAME' ] = 'srel_' + fftype
	d['CMD'] = 'python srel_localdensity.py %s %s %s' % (lammpstraj, fftype, 'c25_'+fftype)
	d['hasMAIL'] = ''
	sJob = 'c25_' + fftype + '_srel' + '.sh'
	sJob = os.path.join(os.path.abspath(target_dir), sJob)
	file(sJob, 'w').write(sJobIn % d)
	os.system('chmod 777 ' + sJob)
	
os.system('cp %s %s' % ('srel_localdensity.py', target_dir))


