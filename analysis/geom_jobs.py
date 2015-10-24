#!/usr/bin/env python

import os
chain_lengths = [25]
sJobIn = file(os.path.expanduser('~/job.template')).read()

lines = open('trajdata.txt', 'r').readlines()
d = {'JOBNAME': 'geom_analyse'}
data_dir = os.path.expanduser('~/c25ld/data')
analysis_dir = os.path.expanduser('~/c25ld/data/analysis')
    
   
for clen in chain_lengths:
    ind = [lines.index(line) for line in lines if line.startswith('%d' % clen)]
    start = ind[0]; stop = ind[-1]+1
    
    # create savedir
    savedir = lines[stop].strip()
    savedir = os.path.join(analysis_dir, savedir)
    if not os.path.isdir(savedir):  os.mkdir(savedir)
    
    # extract params
    cmdstring = ''
    for this_line in lines[start:stop]:
        l = this_line.strip().split()
        nmon = l[0]; nwater = l[1]; npoly = 1
        trajtype = l[2]
        trajfile = l[-1]; trajfile = os.path.join(data_dir, trajfile)
        cmdstring += 'python geom_analyse.py %s "%s" %s %s %s %s\n' % (trajfile, trajtype, nmon, npoly, nwater, savedir)
        print cmdstring
    
    d['CMD'] = cmdstring        
    sJob = 'geom_analyse_%d.sh' % clen
    file(sJob, 'w').write(sJobIn % d)
    os.system('chmod 777 %s' % sJob)
    
