#!/usr/bin/env python
import os, sys, pickle
import numpy as np
import sim
from sim.utility import ProgressBar as pb
import pickleTraj

Normalize = True

## Constant parameters
## The cut-offs are taken on the basis of first-shell radii
## of methane-methane rdfs kept in ~/c25ld/data/analysis/methane25_rdf/
rdfcut = {'lj': 5.5, 'wca': 5.5}  

## User input
traj = os.path.abspath(sys.argv[1])
fftype = sys.argv[2]
nmethanes = int(sys.argv[3])
prefix = sys.argv[4]
savedir = sys.argv[5]

## check if computation needs to be done
pickleName = os.path.join(savedir, prefix + '.pickle')
if os.path.isfile(pickleName):
    quit()

## pickle trajectory
trj = pickleTraj(traj)

## frame-selector
start = 0; stop = len(trj); freq = 10
frame_range = range(start,stop,freq)
nframes = len(frame_range)

## Initialise trajectory reader and read parameters
boxl = trj.FrameData['BoxL']

## make bins
bin_centers = range(0, nmethanes+1)
bin_val = np.zeros(nmethanes+1, np.float64)

## frame stepping
frameStatus = pb(Text = 'Reading frame by frame', Steps = nframes)
for frame in frame_range:
	## use clustering routine kept in sim
	Pos = trj[frame]
	clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = boxl, Cutoff = rdfcut[fftype])
	bin_val += np.array(clustdist)
	frameStatus.Update(frame/freq)

## average over all the frames
bin_val /= nframes

## normalize data
if Normalize:
	bin_val /= np.sum(bin_val)

## pickle data
pickle.dump((bin_centers, bin_val), open(pickleName, 'w'))
