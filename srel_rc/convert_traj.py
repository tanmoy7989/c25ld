#!/usr/bin/env python

import os, pickle, sys
import numpy as np
import sim
import pickleTraj

LammpsTrajIn = os.path.abspath(sys.argv[1])
LammpsTrajOut = os.path.abspath(sys.argv[2])

Trj = pickleTraj(LammpsTrajIn)
BoxL = Trj.FrameData['BoxL']
monomer_atomname = '3'

mappedatoms = []
for (i,AtomType) in enumerate(Trj.AtomTypes):
    if AtomType == int(monomer_atomname):
        mappedatoms.append(i)
    
Map = sim.atommap.PosMap()
for (i,j) in enumerate(mappedatoms):
    Map += [sim.atommap.AtomMap(Atoms1 = j, Atom2 = i, Atom2Name = monomer_atomname)]

AtomNames = [monomer_atomname] * len(mappedatoms)
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomNames)
sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrajOut, Verbose = True)
pickleTraj(LammpsTrajOut)
