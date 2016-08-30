#!/usr/bin/env python
import sim
import matplotlib.pyplot as plt
import numpy as np

data = eval(file('sim.data').read())

# build system in sim
print 'Making system...'
NMol = 100
atomtype = sim.chem.AtomType('M', Mass = 16.0427, Charge = 0.0)
moltype = sim.chem.MolType('Methane', [atomtype])
world = sim.chem.World([moltype], Dim = 3, Units = sim.units.AtomicUnits)
Sys = sim.system.System(world, Name = 'matchEne')
for i in range(NMol): Sys += moltype.New()

p_wca = sim.potential.LJ(Sys, Filter = sim.atomselect.Pairs, Sigma = 3.73, Epsilon = 0.13986, Cut = 4.187, Label = 'P_WCA')
p_ld = sim.potential.LocalDensity(Sys, Filter = sim.atomselect.PolyFilter([atomtype, atomtype], Ordered = True), 
								  Cut = 7.8, LowerCut = 6.6, Knots = data['LDKnots'], RhoMin = 0, RhoMax = 5., Label = 'P_LD')

Sys.ForceField.extend([p_wca, p_ld])
for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
Sys.Measures.PEnergy.SetupHist(-500., 500., 300)
Sys.Load()

Sys.TempSet = 298.0
Sys.BoxL = [35.] * 3
#Sys.Arrays.Pos = data['InitPos']
sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = 298.0)

Int = Sys.Int
for Method in Int.Methods:
	if hasattr(Method, 'TimeStep'): Method.TimeStep = 3.e-3
Sys.Measures.VerboseOutput(StepFreq = 100)

# export to lammps
sim.export.lammps.LammpsExec = 'lmp_tsanyal'
Trj, TrjFile = sim.export.lammps.MakeLammpsTraj(Sys = Sys, Prefix = Sys.Name, NStepsMin = 1000, NStepsEquil = 10000, NStepsProd = 10000, WriteFreq = 100)

# match energies
print 'Computing energies...'
sim_pe = []
lmp_pe = []
diff = []
fracdiff = []
for (i,Pos) in enumerate(Trj):
	Sys.Arrays.Pos = Pos
	Sys.ForceField.Eval()
	sim_pe.append(Sys.PEnergy)
	lmp_pe.append(Trj.FrameData['PEnergy'])
	diff.append(np.abs(sim_pe[-1] - lmp_pe[-1]))
	fracdiff.append(diff[-1]/sim_pe[-1])

N = len(sim_pe)
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1) ; ax2 = fig.add_subplot(2,1,2)
ax1.plot(range(N), sim_pe, linestyle = 'solid', lw = 2, label = 'sim')
ax1.plot(range(N), lmp_pe, linestyle = 'none', marker = 'o', markersize = 5, label = 'lmp')
ax1.set_xlabel('timestep') ; ax1.set_ylabel('total pe')

ax2.plot(range(N), diff)
#ax2.set_ylim([-1.e-6, 1.e-4])
ax2.set_xlabel('timestep') ; ax2.set_ylabel('relative error')

plt.show()






