#!/usr/bin/env python

import sim
import sys
import pickleTraj

LammpsTrjIn = sys.argv[1]
LammpsTrjOut = sys.argv[2]
N_mol = 500
Trj = pickleTraj(LammpsTrjIn, Verbose = False)
BoxL = Trj.FrameData['BoxL']
Map = sim.atommap.PosMap()
for i in range(0, N_mol): Map += [sim.atommap.AtomMap(Atoms1 = i*3, Atom2= i)]
AtomTypes = [1]*N_mol
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes)
#sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrjOut, Verbose = True)
pickleTraj(LammpsTrjOut, Verbose = True)       
