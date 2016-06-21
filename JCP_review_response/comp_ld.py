#!/usr/bin/env python

import os, sys, pickle
import numpy as np

sys.path.append(os.path.expanduser('~/selLDCut'))
import structcorr as ld

ld.calcErrorBar = False
ld.Normalize = True
ld.LDCut = 7.8
ld.LDDelta = 1.2
ld.AtomNames2Types = False
ld.MeasureFreq = 10
ld.Nbins = 100

Nmon = [25, 40, -1]
TrajDir = os.path.expanduser('~/c25ld/data/lammpstraj')

for n in Nmon:
	if n == 25: ld.LammpsTraj = os.path.join(TrajDir, 'c25/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz')
	if n == 40: ld.LammpsTraj = os.path.join(TrajDir, 'polymer_transferability_runs/c40/unconstrained_wca/c40_unbiased_mapped.lammpstrj.gz')
	if n == -1: ld.LammpsTraj = os.path.expanduser('~/c25ld/newcutoff/modtraj/c40/c40_wca_SPLD.lammpstrj.gz')
	if n == -1: ld.Prefix = 'c40_SPLD'
	else: ld.Prefix = 'c%d' % n
	ld.genFileNames()
	ld.makeLD(0,0)