#!usr/bin/env python
import os, sys
import topology as topo

sys.dont_write_bytecode = True

template = "./templates/lammpsdata_template"
start  = 5
atomtype_dict = {"H":1, "O":2, "P":3}
charge_dict = {"H": 0.4238, "O":-0.8476, "P":0.0}


def getData(PDBfile):
	with open(PDBfile,"r") as of:
		s = of.readlines()
	of.close()
	return s


def writeLammpsData(s,fname, BoxPos):
	atomstring = s[0]
	bondstring = s[1]
	anglestring = s[2]

	Natoms = len(atomstring.strip().split("\n"))
	Nbonds = len(bondstring.strip().split("\n"))
	Nangles = len(anglestring.strip().split("\n"))

	#print "Writing to file..........."

	with open(template, "r") as f1:
		src = f1.readlines()
	f2 = open(fname,"w")
	trg = ""
	for line in src:
		if line.startswith("Natoms"): line = line.replace("Natoms", str(Natoms))
		if line.startswith("Nbonds"): line = line.replace("Nbonds", str(Nbonds))
		if line.startswith("Nangles"): line = line.replace("Nangles", str(Nangles))
		if line.startswith("atom_id"): line = line.replace(line, atomstring)
		if line.startswith("bond_id"): line = line.replace(line,bondstring)
		if line.startswith("angle_id"): line = line.replace(line, anglestring)

		#set co-ordinates of bounding box
		if line.find("xhi")>1: line = "%g %g xlo xhi\n" %(BoxPos[0],BoxPos[3])
		if line.find("yhi")>1: line = "%g %g ylo yhi\n" %(BoxPos[1],BoxPos[4])
		if line.find("zhi")>1: line = "%g %g zlo zhi\n" %(BoxPos[2],BoxPos[5])

		f2.write(line)

	f1.close()
	f2.close()


def getmol_id(pdbline):
	if pdbline.split()[3]=="HOH":
		if (pdbline.split()[4]).isalpha():
			mol_id = int(pdbline.split()[5])
		else:
			mol_id = pdbline.split()[4]
			mol_id = mol_id[1:]
			mol_id = int(mol_id)

	if pdbline.split()[3]=="pol":
		mol_id = int(pdbline.split()[5])



	return mol_id


def makeAtoms(pdbstring):
	print "Generating atomlist"
	s = ""
	prev_id = getmol_id(pdbstring[start])
	mol_id = 1;
	for j,line in enumerate(pdbstring[start:len(pdbstring)-1]):
		z = float(line.split()[-1])
		y = float(line.split()[-2])
		x = float(line.split()[-3])
		atom_id = j+1
		curr_id = getmol_id(line)
		mol_id += (not curr_id==prev_id)
		prev_id = curr_id
		atom_type = [atomtype_dict[key] for key in atomtype_dict.keys() if key==line.split()[2]][0]
		atom_type = int(atom_type)
		q = charge_dict[line.split()[2]]
		s += "%d %d %d %f %f %f %f\n" %(atom_id, mol_id, atom_type, q, x, y, z)
	return s



def makeBonds(atomstring):
	print "Generating Bonds"
	s = ""
	nbond = 0
	atomstring = atomstring.strip()
	for i, line in enumerate(atomstring.split("\n")):
		thisatom = line
		for line_nxt in atomstring.split("\n")[i+1:]:
			nextatom = line_nxt
			if thisatom.split()[1] == nextatom.split()[1]:
				(hasBond,btype) = topo.isBonded(thisatom,nextatom)
				if hasBond:
					nbond += 1
					s += "%d %d %d %d\n" %(nbond, btype, int(thisatom.split()[0]), int(nextatom.split()[0]))
	return s



def makeAngles(bondstring):
	print "Generating Angles"
	s = ""
	nangle = 0
	bondstring = bondstring.strip()
	for i, line in enumerate(bondstring.split("\n")):
		thisbond = line
		for line_nxt in bondstring.split("\n")[i+1:]:
			nextbond = line_nxt
			(hasAngle,atype,angle) = topo.isAngle(thisbond,nextbond)
			if hasAngle:
				nangle += 1
				flankleft = angle[0]; center = angle[1]; flankright = angle[2]
				s += "%d %d %d %d %d\n" % (nangle, atype, flankleft, center, flankright)

	return s




