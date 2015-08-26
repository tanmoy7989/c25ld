#!/usr/bin/env python

TempSet = d['TempSet']

RelaxSteps = 100000
    atomtype = sim.chem.AtomType("P", Mass = poly_dict["Mass_mon"], Color = (1,0,0), Radius = 0.5)
    Sys.BoxL[:] = BoxL
    if fftype == 'lj':
    	Sys.ForceField.extend([Plj,Pspline, Pbond, Pangle, Plocaldensity])
    	Sys.ForceField.SetSplineTreatment( NonbondEneSlope = 50., BondEneSlope = 10., AngleEneSlope = 60.)
    	Sys.ForceField.extend([Plj, Pbond, Pangle])
    ## For constrained simulations the biasing potential is added *before* compilation
    ## So, for constrained (umbrella) simulations, return the System obejct at this point
    Int = Sys.Int
    if isUmbrella:
    	return Sys
    	
    ## Load and compile the system
    Sys.Load()
    sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
    for Method in Int.Methods:
    Int.Method = Int.Methods.VVQuench