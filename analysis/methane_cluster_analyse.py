#!/usr/bin/env python
import os, sys, pickle
import numpy as np
import sim
from sim.utility import ProgressBar as pb


doMinImage = True
Normalize = False

## Constant parameters
## The cut-offs are taken on the basis of first-shell radii
## of methane-methane rdfs kept in ~/c25ld/data/analysis/methane
N_mon = 25
N_water = 1700
cut = {'lj': 5.5, 'wca': 5.5}  


## User input
trajname = os.path.abspath(sys.argv[1])
if sys.argv[2]:	start = int(sys.argv[2])
else: start = 0
if sys.argv[3]:	stop = int(sys.argv[3])
else: stop = -1
if sys.argv[4]: freq = int(sys.argv[4])
else: freq = 1
fftype = sys.argv[5]
prefix = sys.argv[6]

## find methane atomtype
cgtype = prefix.split('_')[1]
if cgtype == 'AA':
	methane_atomtype = 3
elif cgtype == 'CG':
	methane_atomtype = 0

## pickle trajectory
if os.path.isfile(trajname + '.pickle'):
	print 'Loading trajectory from pickle...'
	trj = pickle.load(open(trajname + '.pickle', 'r'))
else:
	print 'Pickling trajectory...'
	trj = sim.traj.lammps.Lammps(trajname)
	pickle.dump(trj, open(trajname + '.pickle', 'w'))

## frame-selector
init = trj[0]
if stop <0: stop = len(trj)
frame_range = range(start,stop,freq)
nframes = len(frame_range)


## Initialise trajectory reader and read parameters
boxl = trj.FrameData['BoxL']
ind = np.where(trj.AtomTypes == methane_atomtype)[0]
cutsq = cut[fftype]**2.

## make bins
bin_centers = range(N_mon+1)
bin_val = np.zeros(N_mon+1, np.float64)

## frame stepping
frameStatus = pb(Text = 'Reading frame by frame', Steps = nframes)

for frame in frame_range:
	
	## use clustering routine kept in sim
	Pos = trj[frame][ind]
	clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = trj.FrameData['BoxL'], Cutoff = cut[fftype])
	bin_val += np.array(clustdist)
	
	frameStatus.Update(frame/freq)

## average over all the frames
bin_val /= nframes

## normalize data
if Normalize:
	bin_val /= np.sum(bin_val)

## pickle data
pickleName = prefix + '.pickle'
pickle.dump((bin_centers, bin_val), open(pickleName, 'w'))
