#!/usr/bin/env python

import os, sys, pickle
import parse_potential as pp

sys.path.append('/home/cask0/home/tsanyal/c25ld/c25_single_site_water')
import cgmodel as cg

cg.Prefix = 'methane25_wca_ssw'
Sys = cg.makeSys()
cg.runSrel(Sys)