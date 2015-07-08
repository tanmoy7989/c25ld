#!/usr/bin/env python

import numpy as np
import sys

sys.dont_write_bytecode = True


moltype_dict = {1: "wat", 2: "wat", 3: "poly"}
blength = {"poly": 1.53, "wat": 1.}
btype_dict = {"wat": 1, "poly": 2}
bondtoangle_type_dict = {1:1, 2:2}
tol = 1e-2

def isBonded(atom1,atom2):

	hasBond = 0
	btype = 0

	atomtype_1 = int(atom1.split()[2])
	atomtype_2 = int(atom2.split()[2])
	moltype = moltype_dict[atomtype_1]
	blen = blength[moltype]
	pos1 = [float(atom1.split()[-3]), float(atom1.split()[-2]), float(atom1.split()[-1])]
	pos2 = [float(atom2.split()[-3]), float(atom2.split()[-2]), float(atom2.split()[-1])]
	pos1 = np.array(pos1)
	pos2 = np.array(pos2)
	this_blen = np.linalg.norm(pos1-pos2)
	if abs(this_blen - blen) <= tol:
		btype = btype_dict[moltype]
		hasBond = 1

	return (hasBond, btype)






def isAngle(bond1,bond2):

	hasAngle = 0
	atype = 0
	angle = (0,0,0)
	center = 0

	atoms = {}
	btype1 = int(bond1.split()[1])
	btype2 = int(bond2.split()[1])
	atoms[int(bond1.split()[-1])] = btype1
	atoms[int(bond1.split()[-2])] = btype1
	atoms[int(bond2.split()[-1])] = btype2
	atoms[int(bond2.split()[-2])] = btype2

	idlist = [int(bond1.split()[-2]),int(bond1.split()[-1]),int(bond2.split()[-2]),int(bond2.split()[-1])]

	for i, atom in enumerate(idlist):
		if idlist.count(atom)>1:
				center = atom
				idlist.pop(i)

	if center:
		hasAngle = 1
		for j, atom in enumerate(idlist):
			if atom==center: idlist.pop(j)

		angle = (idlist[0],center,idlist[1])
		atype = atoms[center]

	return(hasAngle,atype,angle)






