#!/usr/bin/env python

import numpy as np
import sys
import sim
 
## contains all different essential subroutines used

def makeSysData(isPolymer = False, hasLJ = False, N_mon = 25, N_water = 1700):
	'''
	Returns dictionaries containing data about the polymer-water system
	    and about the different cutoffs used in the system
    '''
    # Box size for methane simulations = 36 A
	# Pair cutoff used for AA methane simulations = 10 A
	# Cutoff for CG methane needs to be more but also less than half the box width
	# So, newcutoff = 15 A 
	## added on 08/20/2015
	
	N_water = 0 #always using mapped trajectories
	
	if isPolymer:
	   N_poly = 1
	   N_mon = N_mon
	   SPCut = 7.0
	else:
	   N_poly = N_mon
	   N_mon = 1
	   SPCut = 15.0
	      
	poly_dict =      {"N_mon": N_mon, "N_poly": N_poly , "N_water" : N_water,
             	     "blength"      :   1.5300,     "Mass_mon"  :  16.0427,
              	     "Kbond"        :   400,   	     "r0"       :   1.53,
             	     "Kangle"       :  55.17204,    
             	     "theta0"    : 111.00*(np.pi/180.0),
             	     "LJEpsilon"    :  0.13986,     "LJSigma"   :    3.73}

	## LD cut will be supplied to makesys
	cut_dict= {"SPCut": SPCut, "LJCut" : 10.0, "WCACut" : 4.187}
	
	## flags
	flag_dict = {'isPolymer': isPolymer, 'hasLJ': hasLJ}
				
	ret = {'poly_dict': poly_dict, 'cut_dict': cut_dict, 'flag_dict': flag_dict}
	return ret



def makeSys(LDCut, SysData, fftype = 'wca', BoxL = None, Prefix = 'methane_wca', paramstring = None, Delta = 1.2, sample_trajtype = 'CG', eval_trajtype = 'CG'):
    '''
    Creates a system and assembles the potentials
    Most general case for both polymer and methane
    '''
    
    poly_dict = SysData['poly_dict']
    cut_dict = SysData['cut_dict']
    flag_dict = SysData['flag_dict']
    isPolymer = flag_dict['isPolymer']; hasLJ = flag_dict['hasLJ']
    
    atomtype = sim.chem.AtomType("P", Mass = poly_dict["Mass_mon"], Color = (1,0,0), Radius = 0.5)
    moltype = sim.chem.MolType("M", [atomtype]*poly_dict["N_mon"])
    
    # create bonds if polymer
    if isPolymer:
        for i in range(0, poly_dict["N_mon"]-1): moltype.Bond(i,i+1)
    
    world = sim.chem.World([moltype], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(world, Name = Prefix)
    for i in range(poly_dict["N_poly"]):
        Sys += moltype.New()
    
    if BoxL is None:
        Sys.BoxL[:] = 0.
    else:
        Sys.BoxL[:] = BoxL
    
    # add the forcefields in correct order and based on sample and eval traj type
    Sys.__setattr__('sample_trajtype', sample_trajtype)
    Sys.__setattr__('eval_trajtype', eval_trajtype)
    
    # spline and local-density (always there)
    if isPolymer:
        splineFilter = sim.atomselect.NonbondPairs4
    else:
        splineFilter = sim.atomselect.NonbondPairs
    Pspline = sim.potential.PairSpline(Sys, Filter = splineFilter,
                                       Cut = cut_dict["SPCut"], NKnot = 40, Label = "SP")
    
    LDFilter = sim.atomselect.PolyFilter([atomtype, atomtype], Ordered = True)
    LDRhoMin = 0.0
    LDRhoMax = 25
    LDNKnot = int((LDRhoMax - LDRhoMin)+1)
    #LDNKnot = 50 # doubling LDKnot to see if any change on methane potentials - 02/15/2016
    Plocaldensity = sim.potential.LocalDensity(Sys, Filter = LDFilter,
                                               Cut = LDCut, LowerCut = LDCut - Delta,
                                               NKnot = LDNKnot, RhoMin = LDRhoMin, RhoMax = LDRhoMax,
                                               Label = "LD")
    
    # bond and angle potentials for polymer
    if isPolymer:        Pbond = sim.potential.Bond(Sys, Filter = sim.atomselect.BondPairs,                                   Dist0 = poly_dict["r0"], FConst = poly_dict["Kbond"],                                   Label = "Bond_AA")        Pangle = sim.potential.Angle(Sys, Filter = sim.atomselect.BondTriples,                                     Theta0 = poly_dict["theta0"], FConst = poly_dict["Kangle"],                                     Label = "Angle_AA")
    
    # LJ/WCA potentials incase LD-only simulations were done
    if hasLJ or Sys.eval_trajtype == 'AA':
        if fftype == 'lj':            Plj_label = 'LJ_AA'            Plj_cut = cut_dict['LJCut']        else:            Plj_label = 'WCA_AA'            Plj_cut = cut_dict['WCACut']        Plj = sim.potential.LJ(Sys, Filter = sim.atomselect.NonbondPairs4, Cut = Plj_cut,                           Sigma = poly_dict["LJSigma"], Epsilon = poly_dict["LJEpsilon"],                           Shift = True, Label = Plj_label)
    
    
    # add potentials in the correct sequence
    if Sys.eval_trajtype == 'AA':
        Sys.ForceField.extend([Plj])
    elif Sys.eval_trajtype == 'CG':
        if hasLJ: Sys.ForceField.extend([Plj])
        Sys.ForceField.extend([Pspline])
        if isPolymer: Sys.ForceField.extend([Pbond, Pangle])
        Sys.ForceField.extend([Plocaldensity])
        # feed in forcefield if available
        if not paramstring is None: Sys.ForceField.SetParamString(paramstring)
    
    # set up the integrator and time stepping 
    Int = Sys.Int
    timestep = 3.e-3
    for Method in Int.Methods:        if hasattr(Method, "TimeStep"): Method.TimeStep = timestep    for P in Sys.ForceField:    P.Arg.SetupHist(NBin = 10000, ReportNBin = 100) 
  
    Sys.Load()
    
    return Sys
   

def getLDSpline(Sys):
    '''
    extracts spline coefficients for the LD potential
    '''
    # extract the local density potentials preserving correct sequence
    # LD potential is always last in the list so this is easy to do
    Plocaldensity = Sys.ForceField[-1]
    
    LDNKnot = Plocaldensity.NKnot
    SPDist = Plocaldensity.SPDist
    SPC0 = Plocaldensity.SPC0
    SPC1 = Plocaldensity.SPC1
    SPC2 = Plocaldensity.SPC2
    SPC3 = Plocaldensity.SPC3
    
    return (LDNKnot, SPDist, SPC0, SPC1, SPC2, SPC3)
           


def FEP(BE11, BE22, BE12, BE21, Dir = 1, Verbose = True):
    '''
    calculates forward or avg of forward and reverse 
    free energy differences using simple free energy
    perturbation
    '''
    
    DeltaBE = BE11 - BE12
    FE1 = -np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) - DeltaBE.max()
    DeltaBE = BE22 - BE21
    FE2 = np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) + DeltaBE.max()

    if Dir == 1:
        FE = FE1 
    elif Dir == 2:
        FE = 0.5*(FE1 + FE2)
        
    if Verbose:
        print "\t---> %dD Free Energy Perturbation gives: %f" % (Dir, FE)
        
    return FE



def BennettFE(BE11, BE21, BE22, BE12, Verbose = False, 
              Tol = 1.e-10, MaxIter = 10000):
    '''
    Returns -log(Q2/Q1) = beta2 F2 - beta1 F1 using the Bennett method.  
    BEij is the list of beta*E values for system i evaluated in system j.
    '''
    
    fd = lambda x,w : 1./(np.exp(w) + np.exp(x))
    
    
    #check array sizes
    N1 = len(BE11)
    if len(BE12) != len(BE11):
        raise ValueError("BE11 and BE12 arrays are not the same length.")
    N2 = len(BE21)
    if len(BE22) != len(BE21):
        raise ValueError("BE21 and BE22 arrays are not the same length.")        
    if not (np.all(np.isfinite(BE11)) and np.all(np.isfinite(BE12)) 
            and np.all(np.isfinite(BE21)) and np.all(np.isfinite(BE22))):
        if DEBUG: 
            print "BE11 BE12 BE22 BE21"
            for i in xrange(n):
                print i, BE11[i], BE12[i], BE22[i], BE21[i]
        raise ValueError("Found non-finite value in BE11, BE12, BE22, or BE21")
    #fe perturbation for initial guesses
    DeltaBE = BE11 - BE12
    FE1 = -np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) - DeltaBE.max()
    DeltaBE = BE22 - BE21
    FE2 = np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) + DeltaBE.max()
    if Verbose:
        print "Bennett initial guesses: %f, %f" % (FE1, FE2)
    M = np.log(float(N2)/float(N1))
    #setup first points: x is beta2 F2 - beta1 F1, y should be zero
    xo = FE1
    NVals = BE21 - BE22 + xo + M
    DVals = BE12 - BE11 - xo - M
    maxterm = -max(NVals.max(), DVals.max())
    Num = np.sum(fd(NVals+maxterm, maxterm))
    Dem = np.sum(fd(DVals+maxterm, maxterm))
    yo = np.log(Num / Dem)
    x = FE2
    y = 1.
    #loop until tolerance is met
    Iter = 0
    while np.abs(y) > Tol and Iter < MaxIter:
        Iter += 1
        #evaluate y for current x
        NVals = BE21 - BE22 + x + M
        DVals = BE12 - BE11 - x - M
        maxterm = -max(NVals.max(), DVals.max())
        Num = np.sum(fd(NVals+maxterm, maxterm))
        Dem = np.sum(fd(DVals+maxterm, maxterm))
        y = np.log(Num / Dem)
        #predict new x
        xn = (y*xo-yo*x)/(y-yo)
        xo = x
        x = xn
        yo = y
        #print messages
        if Verbose:
            print "Bennett iteration %d: current error is %.3e" % (Iter, np.abs(y))
    #now compute the estimated error
    FE = xo
    FEerr = BennettErr(BE11, BE21, BE22, BE12, FE)
    if Verbose:
        print "Bennett final free energy: %f +- %f" % (FE, FEerr)
    return FE, FEerr
    
        
def BennettErr(BE11, BE21, BE22, BE12, FE):
    '''
    Computes the error in the Bennett calculation, sqrt(var(ln(q2/q1))).
    Based on Shirts et al, PRL 91, 140601 (2003)
    '''
    #compute all "work" measurements
    W = np.concatenate((BE12 - BE11, BE22 - BE21))
    #add free energy (ln (Q2/Q1))
    W = W + FE
    n = len(W)
    c = np.max(np.abs(W))
    terms = 1. / (2.*np.exp(-c) + np.exp(W-c) + np.exp(-W-c))
    err = (np.exp(c) / np.mean(terms) - 4.) / n
    err = np.sqrt(err)
    return err    


def calcCoeff(LDCut, Delta = 1.2):
	'''
	Calculates Local Density Coefficients 
	based on given UpperCut. It is assumed that
	InnerCut = Cut - Delta 
	'''
	InnerCut = LDCut - Delta
	LDCutsq = LDCut*LDCut
	InnerCutsq = InnerCut * InnerCut
	eta = InnerCut/LDCut
	etasq = eta * eta
	denom = (1-etasq)
	denom3 = (1-etasq)**3.
	denom4 = denom*denom3
	detadrc = -InnerCut/LDCutsq
	
	c0 = (1-3*etasq)/denom3
	c2 = (1/LDCutsq) * (6*etasq)/denom3
	c4 = - (1/(LDCutsq*LDCutsq)) * 3*(1 + etasq)/denom3
	c6 = (1/(LDCutsq*LDCutsq*LDCutsq)) * (2/denom3)
	coeff = np.array([c0,c2,c4,c6])
	
	dc0 = (- 12*etasq*eta/denom4) * detadrc
	dc2 =  (12*eta*(1+2*etasq)/denom4) * detadrc/LDCutsq - (2/LDCut)*c2	
	dc4 = - (12*eta*(2+etasq)/denom4) * detadrc/(LDCutsq*LDCutsq) - (4/LDCut)*c4
	dc6 = (12*eta/denom4) * detadrc/(LDCutsq*LDCutsq*LDCutsq) - (6/LDCut)*c6
	dcoeff = np.array([dc0,dc2,dc4,dc6])
	
	return (coeff, dcoeff)
	


