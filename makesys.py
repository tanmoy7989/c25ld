#!/usr/bin/env python
import numpy as npimport osimport sysimport pickleimport simfrom chem_data import getModelDatad = getModelData()poly_dict = d['poly_dict']cut_dict = d['cut_dict']
TempSet = d['TempSet']

RelaxSteps = 100000def MakeSys(fftype, BoxL, N_mon = 25, N_poly = 1, Prefix = 'FullSys',           isTest = False, isUmbrella = False):    """    Makes a Sys object based on LJ, Spline, Bond, Angle and Local Density potentials.    If it's a test run then only LJ, Bond and Angle potentials are generated.    """    np.random.seed(12345)    print 'Making system object...'    ## Making the world and the system
    poly_dict['N_mon'] = N_mon
    poly_dict['N_poly'] = N_poly
    atomtype = sim.chem.AtomType("P", Mass = poly_dict["Mass_mon"], Color = (1,0,0), Radius = 0.5)    moltype = sim.chem.MolType("M", [atomtype] * poly_dict["N_mon"])    for i in range(0, poly_dict["N_mon"]-1):    moltype.Bond(i,i+1)    world = sim.chem.World([moltype], Dim = 3, Units = sim.units.AtomicUnits)    Sys = sim.system.System(world, Name = Prefix)    for i in range(poly_dict["N_poly"]): Sys += moltype.New()    ## Constructing the simulation box
    Sys.BoxL[:] = BoxL    ## Desiging the potentials
    if fftype == 'lj':        Plj_label = 'LJ_AA'        Plj_cut = cut_dict['LJCut']    else:        Plj_label = 'WCA_AA'        Plj_cut = cut_dict['WCACut']    Plj = sim.potential.LJ(Sys, Filter = sim.atomselect.NonbondPairs4, Cut = Plj_cut,                            Sigma = poly_dict["LJSigma"], Epsilon = poly_dict["LJEpsilon"],                            Shift = True, Label = Plj_label)    Pspline = sim.potential.PairSpline(Sys, Filter = sim.atomselect.NonbondPairs4,                               Cut = cut_dict["SPCut"], NKnot = 40, Label = "SP")    Pbond = sim.potential.Bond(Sys, Filter = sim.atomselect.BondPairs,                       Dist0 = poly_dict["r0"], FConst = poly_dict["Kbond"],                       Label = "Bond_AA")    Pangle = sim.potential.Angle(Sys, Filter = sim.atomselect.BondTriples,                         Theta0 = poly_dict["theta0"], FConst = poly_dict["Kangle"],                         Label = "Angle_AA")    LDFilter = sim.atomselect.PolyFilter([atomtype, atomtype], Ordered = True)    LDRhoMin = 0.0    LDRhoMax = max(25, poly_dict["N_mon"]) if poly_dict["N_mon"] <=25 else min(25, poly_dict["N_mon"])     LDNKnot = int((LDRhoMax - LDRhoMin) + 1)    Plocaldensity = sim.potential.LocalDensity2(Sys, Filter = LDFilter,                                               Cut = cut_dict["LDCut"], LowerCut = cut_dict["LDLowerCut"],                                               NKnot = LDNKnot, RhoMin = LDRhoMin, RhoMax = LDRhoMax,                                               Label = "LD")    if not isTest:
    	Sys.ForceField.extend([Plj,Pspline, Pbond, Pangle, Plocaldensity])
    	Sys.ForceField.SetSplineTreatment( NonbondEneSlope = 50., BondEneSlope = 10., AngleEneSlope = 60.)    else:
    	Sys.ForceField.extend([Plj, Pbond, Pangle])    ## Compile and load the system only if not a constrained simulation
    ## For constrained simulations the biasing potential is added *before* compilation
    ## So, for constrained (umbrella) simulations, return the System obejct at this point
    Int = Sys.Int    Int.Method = Int.Methods.VVIntegrate    for P in Sys.ForceField:    P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)    
    if isUmbrella:
    	return Sys
    	
    ## Load and compile the system
    Sys.Load()    ## Initialization
    sim.system.init.positions.CubicLattice(Sys, Random = 0.1)    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)        #setup time step
    for Method in Int.Methods:        if hasattr(Method, "TimeStep"): Method.TimeStep = 3.e-3    #run for a bit to minimize and remove clashes
    Int.Method = Int.Methods.VVQuench    Int.Run(RelaxSteps, "Minimizing")    return Sys
