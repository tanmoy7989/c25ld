#!/usr/bin/env python

import os, pickle
import numpy as np
import sim
import pickleTraj

LammpsTrajIn = os.path.expanduser('~/c25ld/data/lammpstraj/methane/nvt/methane_wca_unbiased.lammpstrj.gz')
LammpsTrajOut = os.path.expanduser('~/c25ld/data/lammpstraj/methane/nvt/methane_wca_unbiased_mapped.lammpstrj.gz')
trj = pickleTraj(LammpsTrajIn)
boxlen = trj.FrameData['BoxL']
mappedatoms = []
for (i,AtomType) in enumerate(trj.AtomTypes):
    if AtomType == 3:
        mappedatoms.append(i)
    
Map = sim.atommap.PosMap()
for (i,j) in enumerate(mappedatoms):
    Map += [sim.atommap.AtomMap(Atoms1 = j, Atom2 = i, Atom2Name = '3')]


AtomNames = [3] * len(mappedatoms)
MappedTrj= sim.traj.Mapped(trj, Map, AtomNames = AtomNames)
sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrajOut, Verbose = True)
pickleTraj(LammpsTrajOut)
