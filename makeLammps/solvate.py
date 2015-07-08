#!/usr/bin/env python

import os
import sys
import numpy as np
template_dir = "./templates"

sys.dont_write_bytecode = True

def createPolymerPos(N_mon, N_poly):
	l = 1.53
	theta = (np.pi / 180.0) * 111.0 / 2
	pos = []
	pos0 = np.random.random((1,3))
	# patch: only 1 molecule creation is required
	# Over-ride N_poly
	N_poly = 1
	for n_mol in range(N_poly):
		pos0 = pos0 + 5*np.random.random((1,3));
		mol_id = n_mol + 1
		for n_atom in range(N_mon):
			x = pos0[0][0]
			y = pos0[0][1] + n_atom * l * np.sin(theta)
			z = pos0[0][2] + (not np.mod(n_atom, 2) == 0) * l * np.cos(theta)
			atom_dict = {"mol_id": mol_id, "x":x, "y":y, "z":z}
			pos.append(atom_dict)


	return pos



def createPolymerPDB(pos):
	s= ""
	for i,r in enumerate(pos):
		atom_id = i+1
		atom_name = "P"
		mol_name = "pol"
		mol_id = r["mol_id"]
		q = 0.0
        	s+="%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f\n" %(
        	"HETATM", atom_id, atom_name, "", mol_name, "", mol_id, "",
        	r["x"], r["y"], r["z"])

	with open('polymer.pdb','w') as of:
    		of.write(s)


#fixed box solvation routine
def solvate_fbox(polyPDB, waterPDB, N_poly,N_water,BoxL):
	f1 = open(template_dir+"/packmol_template", "r")
	f2 = open("packmol_input.inp", "w")
	for line in f1:
		if line.endswith("waterfile\n"): line = line.replace("waterfile", waterPDB)
		if line.endswith("polyfile\n"): line = line.replace("polyfile", polyPDB)
		if line.endswith("N_water\n"): line = line.replace("N_water", str(N_water))
		if line.endswith("N_poly\n"): line = line.replace("N_poly", str(N_poly))
		if line.endswith("Len\n"): line = line.replace("Len", str(BoxL/2.))
		f2.write(line)
	f1.close()
	f2.close()
	os.system("packmol<packmol_input.inp")



#fixed shell solvation routine (uses automatic solvator tcl script)
def solvate_fshell(polyPDB,solv_shell):

	os.system("solvate.tcl "+ polyPDB + " -shell " + str(solv_shell) + " -noions -density 1.00 -o solvated_poly.pdb")

	# parse the waterpdb file to match rest of code
	f1 = open("WATER.pdb", "r")
	f2 = open("water.pdb", "w")
	for line in f1:
		if line.split()[2] == "H1":
			line=line.replace("H1" ,"H")
			line =  line.replace("TIP3", "HOH")
		if line.split()[2] == "H2":
			line=line.replace("H2" ,"H")
			line =  line.replace("TIP3", "HOH")
		if line.split()[2] == "OH2":
			line=line.replace("OH2" ,"O")
			line =  line.replace("TIP3", " HOH")

		f2.write(line)
	f1.close()
	f2.close()
	os.remove("WATER.pdb"); os.system("mv water.pdb WATER.pdb")
	os.system("packmol<packmol_input.inp")



def getBoxPos(packmol_file, solv_type):
	# get box co-ordinates and other number data for lammps input file

	lines = open(packmol_file, "r").readlines()
	solvtype_dict = {"fixed_box":["structure " + template_dir +"/water.pdb\n", 3],
		     "fixed_shell":["structure WATER.pdb\n", 2]}

	if lines.__contains__(solvtype_dict[solv_type][0]): windex = lines.index(solvtype_dict[solv_type][0])
	gap = solvtype_dict[solv_type][1]

	N_water  = int(lines[windex+1].split()[-1])	#get number of water
	pos = lines[windex+gap].split()[2:]		#get box pos
	pos = [float(k) for k in pos]

	return(N_water,pos)



def correctBoxPos(boxpos,N_mon):
	#correct box pos to avoid overlaps between images
	l = 1.53
	maxL = (N_mon*l)
	x0 = boxpos[0];x1 = boxpos[3]
	y0 = boxpos[1];y1 = boxpos[4]
	z0 = boxpos[2];z1 = boxpos[5]
	pos = [(x0,x1), (y0,y1), (z0,z1)]

	#adding tolerance of 2 Angstroms to boxlength is standard PBC requirement suggested in PackMol
	correctpos = [pos[0][0]-1., pos[1][0]-1., pos[2][0]-1., pos[0][1]+1., pos[1][1]+1., pos[2][1]+1.]

	print "L/2 >= ", maxL
	print "initial pos ", boxpos
	print "corrected pos ", correctpos
	return correctpos


def cleanup():
	redundant_list = ["WATER.pdb", "CHLORIDE.pdb", "SODIUM.pdb", "polymer.pdb", "packmol_input.inp"]
	for file in redundant_list:
		if os.path.isfile(file): os.remove(file)





####test#####



