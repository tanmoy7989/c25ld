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

def getMonIndex(posIndex, indtype = 'relative'):
	if indtype == 'relative':
		return posIndex - (3*N_water)
	elif indtype == 'absolute':
		return posIndex + (3*N_water)


def doMinImage(pos, boxl):	
	for i in range(1, N_mon):
		pos[i] = sim.geom.Reimage(pos[i], pos[i-1], boxl)
	return pos

		
def isNeigh(pos1, pos2, cutsq):
	rsq = np.sum(pos1**2. - pos2**2.)
	return rsq<=cutsq
	


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
bin_centers = range(1,N_mon+1)
bin_val = np.zeros([N_mon], np.float64)

## frame stepping
frameStatus = pb(Text = 'Reading frame by frame', Steps = nframes)

for frame in frame_range:
	## initialize methane_list and positions
	methane_list = range(N_mon)
	pos = trj[frame][ind]
	
	## minimum image
	if doMinImage:
		pos = doMinImage(pos, boxl)
	
	## populate hashmap
	for i in methane_list:
		# initialize on-the-fly hashmap
		clust_size = 1
		clust_members = [i]
		for j in methane_list:
			if j <= i: continue
			## check neighborhood 
			if isNeigh(pos[i], pos[j], cutsq):
				clust_size += 1
				clust_members.append(j)
		
		## update methane_list and bin_count
		[methane_list.remove(x) for x in clust_members]
		bin_val[clust_size-1] += 1.
		
	frameStatus.Update(frame/freq)

## average over all the frames
bin_val /= nframes

## normalize data
if Normalize:
	bin_val /= np.sum(bin_val)

## pickle data
pickleName = prefix + '.pickle'
pickle.dump((bin_centers, bin_val), open(pickleName, 'w'))
