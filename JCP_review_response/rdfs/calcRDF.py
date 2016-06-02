#!/usr/bin/env python

import os, sys
import numpy as np

sys.path.append(os.path.expanduser('~/selLDCut'))
from selLDCut import structcorr as rdf

LammpsTraj = sys.argv[1]
Prefix = sys.argv[2]

rdf.AtomNames2Types = False
rdf.calcErrorBar = False
rdf.Nbins = 50
rdf.TrjIter = [0, -1, 10]

rdf.LammpsTraj = LammpsTraj
rdf.Prefix = Prefix
rdf.genFileNames()
rdf.makeRDF([0], [0])
