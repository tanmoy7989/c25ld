#!/usr/bin/env python

import sim
import sys
import pickleTraj

LammpsTrjIn = sys.argv[1]
LammpsTrjOut = sys.argv[2]
N_mon = int(sys.argv[3])
N_water = int(sys.argv[4])

Trj = pickleTraj(LammpsTrjIn, Verbose = False)
BoxL = Trj.FrameData['BoxL']

Map = sim.atommap.PosMap()
for i in range(0, N_water): Map += [sim.atommap.AtomMap(Atoms1 = i*3+2, Atom2= i)]
for i in range(0, N_mon):   Map += [sim.atommap.AtomMap(Atoms1 = 3*N_water + i, Atom2 = N_water + i)]

AtomTypes = [1]*N_water + [2]*N_mon
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes)
sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrjOut, Verbose = True)       
