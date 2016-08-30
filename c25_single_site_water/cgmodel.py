import os, sys, pickle
from copy import copy
import parse_potential as pp
import pickleTraj

sys.path.append('/home/cask0/home/tsanyal/software')
import mysim
import sim

fixWW = True

LammpsExec = 'lmp_tsanyal'
sim.srel.base.DiffEneFracTol = 0.1
sim.export.lammps.InnerCutoff = 0.05

LammpsTraj = '/home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_single_site_water.lammpstrj.gz'

def makeSys(LDCut = 7.8, NMethane = 25, NWater = 1700):
	print 'Making system'
	atomtypeM = sim.chem.AtomType('M', Mass = 16.0427, Charge = 0.0)
	atomtypeW = sim.chem.AtomType('W', Mass = 18.0, Charge = 0.0)
	moltypeM = sim.chem.MolType('M', [atomtypeM])
	moltypeW = sim.chem.MolType('W', [atomtypeW])
	world = sim.chem.World([moltypeM, moltypeW], Dim = 3, Units = sim.units.AtomicUnits)
	Sys = sim.system.System(world, Name = 'findtimestep')
	for i in range(NWater): Sys += moltypeW.New()
	for i in range(NMethane): Sys += moltypeM.New()

	Filter_MM = sim.atomselect.PolyFilter([atomtypeM, atomtypeM], MinBondOrd = MinBondOrd) 
    Filter_MW = sim.atomselect.PolyFilter([atomtypeM, atomtypeW], MinBondOrd = MinBondOrd)
    Filter_WW = sim.atomselect.PolyFilter([atomtypeW, atomtypeW], MinBondOrd = MinBondOrd)
    Filter_MW_LD = sim.atomselect.PolyFilter([atomtypeM, atomtypeW], Ordered = True)
    Delta = 1.2

    Pspline_MM = sim.potential.PairSpline(Sys, Filter = Filter_MM, Cut = 10.0, NKnot = 30, Label = "SP_MM")
    Pspline_MW = sim.potential.PairSpline(Sys, Filter = Filter_MW, Cut = 10.0, NKnot = 30, Label = "SP_MW")
    Pspline_WW = sim.potential.PairSpline(Sys, Filter = Filter_WW, Cut = 10.0, NKnot = 30, Label = "SP_WW")
    Plocaldensity_MW = sim.potential.LocalDensity(Sys, Filter = LDFilter, Cut = LDCut, LowerCut = LDCut - Delta, NKnot = 30, RhoMin = 0., RhoMax = NMethane, Label = "LD")

    if fixWW:
    	WWSPKnots = eval(file(os.path.expanduser('~/c25ld/c25_single_site_water/water_spce/spce_ff.dat')).read())
    	Pspline_WW.SetParams(Knots = WWSPKnots)

    Sys.ForceField.extend([Pspline_MM, Pspline_MW, Pspline_WW, Plocaldensity_MW])
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 1000, ReportNBin = 100)
    #Sys.ForceField.SetSplineTreatment(NonbondEneSlope = 50., BondEneSlope = 10., AngleEneSlope = 60.)
    for Pspline in [Pspline_MM, Pspline_MW, Pspline_WW]: Pspline.EneSlopeInner = None

    Sys.Load()

	Trj = pickleTraj(LammpsTraj)
	Sys.BoxL = Trj.FrameData['BoxL']
	Sys.TempSet = 298.0
	sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
	sim.system.init.velocities.Canonical(Sys, Temp = 298.0)

	Int = Sys.Int
	for Method in Int.Methods:
		if hasattr(Method, 'TimeStep'): Method.TimeStep = 3.e-3

	return Sys


def runSrel(Sys):
	Prefix = Sys.Name

    Trj = pickleTraj(LammpsTraj, Verbose = True)
    Map = sim.atommap.PosMap()
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)] #1:1 mapping since using a mapped traj
    
    Int = Sys.Int
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = 0.01
    Sys.TempSet = 298.0

    Pspline_MM = Sys.ForceField[0] ; Pspline_MW = Sys.ForceField[1] ; Pspline_WW = Sys.ForceField[2]
    Plocaldensity = Sys.ForceField[-1] 
    PairPotentialList = [Pspline_MM, Pspline_MW] if fixWW else [Pspline_MM, Pspline_MW, Pspline_WW]
    print "Resetting all potentials to base case"
    for Pspline in PairPotentialList: Pspline.SetParam(Knots = 0.)
    Plocaldensity.SetParam(Knots = 0.)
    for p in Sys.ForceField: p.FreezeParam()

    print Sys.ForceField.ParamString()
  
    Opt = sim.srel.OptimizeTrajLammpsClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix, Verbose = True)
    sim.srel.optimizetraj.PlotFmt = 'svg'
    Opt.MinReweightFrac = 0.2
    
    opt_cases = ['SP', 'SPLD']
    for i, case in enumerate(opt_cases):
        Opt.Reset()
        Opt.FilePrefix = Prefix + '_' + case
        print "\n\nOptimizing for the case: ", case
        if case == "SP":
            for Pspline in PairPotentialList: Pspline.UnfreezeParam()
        if case == "SPLD":
            Plocaldensity.UnfreezeParam()
            
        Sys.ForceField.Update()
        Opt.RunConjugateGradient(StepsEquil = 1000000, StepsProd = 2000000, StepsStride = 500) 
