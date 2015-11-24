#!/usr/bin/env python

import numpy as np
import os, sys
import sim
import pickleTraj

# Globals
trans_CGdir = os.path.abspath(sys.argv[1])
trans_AAdir_fmt = sys.argv[2]
c25ffdir = os.path.abspath(sys.argv[3])
c_lens = [10, 20, 30, 40, 50]
sJobIn = file(os.path.expanduser('~/job.template')).read()

fftypes = ['wca']
cgtypes = ['SP', 'SPLD', 'LD']

for fftype in fftypes:
    for c_len in c_lens:
        AAdir = os.path.abspath(os.path.join(trans_AAdir_fmt % c_len, 'unconstrained_%s' % fftype))
        CGdir = os.path.join(trans_CGdir, 'c%d' % c_len)
        if not os.path.isdir(CGdir): os.mkdir(CGdir)
        os.system('cp cg_MD.py %s' % CGdir)
       
        # get BoxL       
        AAtraj = os.path.join(AAdir, 'c%d_unbiased.lammpstrj.gz' % c_len)
        trj = pickleTraj(AAtraj, Verbose = False)
        boxlen = trj.FrameData['BoxL'][0]
        
        # write job files
        for cgtype in cgtypes:
            ffield_file = os.path.join(c25ffdir, 'c25_%s_%s_sum.txt' % (fftype, cgtype))
            Prefix = 'c%d_%s_%s' % (c_len, fftype, cgtype)
            cmdstring = 'python cg_MD.py %d %s %s %s %g' % (c_len, ffield_file, Prefix, fftype, boxlen)
            sJob = os.path.join(CGdir, '%s_MD.sh' % Prefix)
            file(sJob, 'w').write(sJobIn % {'JOBNAME': 'md_srel_rc', 'CMD': cmdstring})
            os.system('chmod 777 %s' % sJob)
            os.system('qsub %s' % sJob)
