import os, sys, pickle
from copy import copy
import parse_potential as pp
import pickleTraj

sys.path.append('/home/cask0/home/tsanyal/software')
import mysim
import sim

sys.path.append('/home/cask0/home/tsanyal/c25ld/c25_single_site_water')
import cgmodel as cg

LammpsTraj = '/home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_single_site_water.lammpstrj.gz'

# build system
NMethane = 25
NWater = 1700
print 'Making system'
atomtypeM = sim.chem.AtomType('M', Mass = 16.0427, Charge = 0.0)
atomtypeW = sim.chem.AtomType('W', Mass = 18.0, Charge = 0.0)
moltypeM = sim.chem.MolType('M', [atomtypeM])
moltypeW = sim.chem.MolType('W', [atomtypeW])
world = sim.chem.World([moltypeM, moltypeW], Dim = 3, Units = sim.units.AtomicUnits)
Sys = sim.system.System(world, Name = 'findtimestep')
for i in range(NMethane): Sys += moltypeM.New()
for i in range(NWater): Sys += moltypeW.New()

# load the original pair WCA for methane-methane and methane-water
# and the spce water forcefield for water-water
P_MM = sim.potential.LJ(Sys, Filter = sim.atomselect.PolyFilter([atomtypeM, atomtypeM]), Sigma = 3.73, Epsilon = 0.13986, Cut = 4.187, Shift = True, Label = 'P_MM')
P_MW = sim.potential.LJ(Sys, Filter = sim.atomselect.PolyFilter([atomtypeM, atomtypeW]), Sigma = 3.73, Epsilon = 0.13986, Cut = 4.187, Shift = True, Label = 'P_MW')
WWSPKnots = eval(file(os.path.expanduser('~/c25ld/c25_single_site_water/water_spce/spce_ff.dat')).read())
P_WW = sim.potential.PairSpline(Sys, Filter = sim.atomselect.PolyFilter([atomtypeW, atomtypeW]), Knots = WWSPKnots, Cut = 10.0, Label = 'P_WW')
Sys.ForceField.extend([P_MM, P_MW, P_WW])

for P in Sys.ForceField: P.Arg.SetupHist(NBin = 1000, ReportNBin = 100)
Sys.Load()

# other system params
Trj = pickleTraj(LammpsTraj)
Sys.BoxL = Trj.FrameData['BoxL']
Sys.TempSet = 298.0
sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = 298.0)

print Sys.ForceField.ParamString()

# find timestep
sim.integrate.velverlet.FindTimeStepLammps(Sys, NSteps = 2000000, GuessTimeStep = 3.e-3)

