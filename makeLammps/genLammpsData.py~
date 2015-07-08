#!/usr/bin/env python

import numpy as np
import sys
import os
import solvate as sol

sys.dont_write_bytecode = True


template_dir =  "./templates"


### User supplied data
sys_temp = sys.argv[1]
n_water = sys.argv[2]
n_mon = sys.argv[3]
n_poly = sys.argv[4]
Prefix = sys.argv[5]
datafilepath = sys.argv[6]
datafilename = os.path.join(datafilepath, Prefix + ".data")


#### create polymer
blength = 1.53 #monomer-monomer bond length
poly_dict= {"N_poly":int(n_poly), "N_mon":int(n_mon), "N_water":int(n_water)}
print "Creating polymer position array..."
pos = sol.createPolymerPos(poly_dict["N_mon"], poly_dict["N_poly"])
sol.createPolymerPDB(pos)


### solvate polymer with fixed box solvation routine
solv_type = "fixed_box"
print "Solvating polymer with %s routine..." % solv_type
BoxL = 2*poly_dict["N_mon"]*blength + 4.00 # Min-image criteria satisfied and extra 4.0 A padding added
sol.solvate_fbox("polymer.pdb", template_dir+"/water.pdb", poly_dict["N_poly"] , poly_dict["N_water"],BoxL)
(N_water, boxpos) = sol.getBoxPos("packmol_input.inp", solv_type)
boxpos = sol.correctBoxPos(boxpos, poly_dict["N_mon"])
sol.cleanup()


#### Call Lammps Input generator module
import makeLammpsData as makeLammps
pdbfile = os.path.join(os.getcwd(), "solvated_poly.pdb")
pdbstring = makeLammps.getData(pdbfile)
atomstring = makeLammps.makeAtoms(pdbstring)
bondstring = makeLammps.makeBonds(atomstring)
anglestring = makeLammps.makeAngles(bondstring)
makeLammps.writeLammpsData((atomstring,bondstring,anglestring), datafilename, boxpos)


#### clean things up 
print "Cleaning up..."
os.remove("solvated_poly.pdb")


