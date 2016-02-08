#!/usr/bin/env python

import numpy as np
import os, pickle, sys
import sim
import pickleTraj
import parse_potential

#import subroutine module
sys.path.append('/home/cask0/home/tsanyal/srel_rc')
import alg

#flags
DEBUG = False
RefreshLogs = True
DelFortCode = True

#fortran source code
fortcode = 'dSreldrc_fort.f90'
fortlib = 'dSreldrc_fortlib'
fcompiler = 'gfortran'
fort_source = '''
SUBROUTINE getdUdrc(Pos, Nmon, BoxL, LDCut, coeff, dcoeff, EneDeriv, &
				    SPDist, SPC1, SPC2, SPC3, NSPDist, NSPC, NKnot)
	IMPLICIT NONE
	
	REAL(8), INTENT(IN), DIMENSION(0:Nmon-1, 0:2) :: Pos
	INTEGER, INTENT(IN) :: Nmon
	REAL(8), INTENT(IN) :: BoxL, LDCut
	REAL(8), INTENT(IN), DIMENSION(0:3) :: coeff, dcoeff
	
	REAL(8), INTENT(OUT) :: EneDeriv
	
	INTEGER, INTENT(IN) :: NSPDist, NSPC, NKnot
	REAL(8), INTENT(IN), DIMENSION(0:NSPDist-1) :: SPDist
	REAL(8), INTENT(IN), DIMENSION(0:NSPC - 1) :: SPC1, SPC2, SPC3
	
	INTEGER :: i,j, spline_i
	REAL(8) :: rsq, invBoxL, LDCutsq, InnerCutsq, delta, phi, dphidrc
	REAL(8) :: spline_x, spline_t
	REAL(8), DIMENSION(0:2) :: Posi, rij
	REAL(8), DIMENSION(0:Nmon-1) :: rho, dfdrho
	
	!zero out rho and dfdrho array
	DO i = 0, Nmon-1
	   rho(i) = 0.0
	   dfdrho(i) = 0.0
	ENDDO
	
	delta = 1.2
	invBoxL = 1.0/BoxL
	LDCutsq = LDCut * LDCut
	InnerCutsq = (LDCut - delta)*(LDCut - delta)
	spline_x = 0.0
	spline_i = 0
	spline_t = 0.0

	! precompute local densities --> rho
	DO i = 0, Nmon - 2
		! precompute central atom for speed
		Posi = Pos(i,:)
		DO j = i+1, Nmon -1			
			! compute pair-distances
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)	
			
			! compute indicator function
			IF (rsq >= LDCutsq) THEN
				phi = 0.0
			ELSE IF (rsq <= InnerCutsq) THEN
				phi = 1.0
			ELSE
				phi = coeff(0) + rsq*(coeff(1) + rsq*(coeff(2) + rsq*coeff(3)))
			ENDIF
			
			rho(i) = rho(i) + phi
			rho(j) = rho(j) + phi
		ENDDO
	ENDDO
	
	! precompute local density derivatives --> dfdrho
	! spline interpolation to find dfdrho(i) at rho(i), taken from python code in sim package by MSS
	DO i = 0, Nmon - 1
	   spline_x = (rho(i) - MINVAL(rho)) * SPDist(1)
	   spline_i = MAX(MIN(INT(spline_x), NKnot - 2), 0)
	   spline_t = MAX(MIN(spline_x - FLOAT(spline_i), 1.0), 0.0)
	   dfdrho(i) = SPDist(1) * (SPC1(spline_i) + spline_t * (2.0*SPC2(spline_i) + spline_t*3.0*SPC3(spline_i)))
	ENDDO
	
	! now compute <dU_CG/drc>
    EneDeriv = 0.0
	DO i = 0, Nmon - 1
	    dphidrc = 0.0
		! precompute central atom for speed
		Posi = Pos(i,:)
		DO j = 0, Nmon - 1
			IF (i==j) THEN
				CYCLE
			ENDIF
			
			! compute pair-distances
			rij = Pos(j,:) - Posi
			rij(:) = rij(:) - BoxL * DNINT(invBoxL * rij(:)) ! minimum-imaging
			rsq = SUM(rij*rij)	
			
			! compute derivative of indicator function
			IF (rsq >= InnerCutsq .AND. rsq <= LDCutsq) THEN
				dphidrc = dphidrc + dcoeff(0) + rsq*(dcoeff(1) + rsq*(dcoeff(2) + rsq*dcoeff(3))) 
			ENDIF	
		ENDDO
		EneDeriv = EneDeriv + dphidrc * dfdrho(i)	
	ENDDO

END SUBROUTINE
'''

#compiling fortran code
def fcompile():
    cmdstring = 'f2py -c -m %s %s --fcompiler=%s' % (fortlib, fortcode, fcompiler)
    file(fortcode, 'w').write(fort_source)
    os.system(cmdstring)
    if DelFortCode:
        for this_file in [fortcode, fortcode.split('.')[0]+'.o']:
            if os.path.isfile(this_file):   os.remove(this_file)

    return fortlib



####### MAIN #######

'''USAGE: 
sample command line call 
python dSreldrc.py lj 0 0 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/methane/nvt/methane_lj_unbiased_mapped.lammpstrj.gz methane_%s_%d_SPLD_sum.txt methane_%s_%d_MD.lammpstrj.gz methane_lj 
'''

#user input
fftype = sys.argv[1]
isPolymer = bool(int(sys.argv[2]))
hasLJ = bool(int(sys.argv[3]))
N_mon = int(sys.argv[4])
N_water = int(sys.argv[5])
LammpsTraj = sys.argv[6] 
ff_fmtstring = sys.argv[7]
md_fmtstring = sys.argv[8]
datalog_prefix = sys.argv[9]

#system parameters
SysData = alg.makeSysData(isPolymer = isPolymer, hasLJ = hasLJ, N_mon = N_mon, N_water = N_water)
poly_dict = SysData['poly_dict']
TempSet = 298.0 ; kB = 0.001987; beta = 1./(kB * TempSet)
Delta = 1.2

#input files
#base_fmt = {'cg_ff': 'methane_%s_%d_SPLD_sum.txt', 'cg_trj': 'methane_%s_%d_MD.lammpstrj.gz'}
base_fmt = {'cg_ff': ff_fmtstring, 'cg_trj': md_fmtstring}

#write output headers
LDCuts = np.arange(4.5, 10.0, 0.1)
datalog = '%s_dSreldrc.dat' % datalog_prefix
if not os.path.isfile(datalog) or RefreshLogs:
    of = open(datalog, 'w')
    of.write('#LDCut \t dSreldrc\n')
    of.close()

#various indexing schemes
TrjIter = {'AA': [0, -1, 10], 'CG': [0, -1, 10]}

# compile fortcode
fcompile()
from dSreldrc_fortlib import getdudrc
Trj_AA = pickleTraj(LammpsTraj)
BoxL = Trj_AA.FrameData['BoxL'][0]

#begin looping over different LD cutoffs
for i, LDCut in enumerate(LDCuts):
	
	print '\n\nLDCut = %g\n' % LDCut
	
	#calculate coefficients
	(coeff, dcoeff) = alg.calcCoeff(LDCut = LDCut, Delta = 1.2)

	#get AA and CG trj
	Trj_CG = pickleTraj(base_fmt['cg_trj'] % (fftype, i))
	
	#get Spline Coeffs for CG forcefield
	ffield_CG = base_fmt['cg_ff'] % (fftype, i)
	paramstring = parse_potential.parseParamString(ffield_CG)
	Sys = alg.makeSys(LDCut = LDCuts[i], Delta = Delta, BoxL = BoxL, fftype = fftype, Prefix = 'ldsplinecoeff',
                      paramstring = paramstring, sample_trajtype = 'CG', eval_trajtype = 'CG', SysData = SysData)

	(nknot, spdist, spc0, spc1, spc2, spc3) = alg.getLDSpline(Sys = Sys)

	for key in TrjIter.keys():
	   if TrjIter[key][1] == -1:
	       if key == 'AA': TrjIter[key][1] = len(Trj_AA)
	       if key == 'CG': TrjIter[key][1] = len(Trj_CG)
	  
	FrameRange_AA = range(TrjIter['AA'][0], TrjIter['AA'][1], TrjIter['AA'][2])
	FrameRange_CG = range(TrjIter['CG'][0], TrjIter['CG'][1], TrjIter['CG'][2])
	NFrames_AA = len(FrameRange_AA); NFrames_CG = len(FrameRange_CG)
	
	#frame stepping
	####all-atom
	print '---> Calculating AA deriviative...\n'
	EneAA_deriv = 0
	pb = sim.utility.ProgressBar(Text = 'Processing AA frames', Steps = NFrames_AA)
	Ind = 0
	for frame in FrameRange_AA:
	   Pos_AA = Trj_AA[frame]
	   x = 0.0
	   x = getdudrc(pos = Pos_AA, nmon = poly_dict['N_mon']*poly_dict['N_poly'],
				    boxl = BoxL, ldcut = LDCut,coeff = coeff, dcoeff = dcoeff,
					nknot = nknot, spdist = spdist, spc1 = spc1, spc2 = spc2, spc3 = spc3)
					
	   if DEBUG: print "AA->", x
	   EneAA_deriv += x							    																			
	   
	   pb.Update(Ind)
	   Ind += 1
		
		
	####coarse-grained
	del x
	print '---> Calculating CG deriviative...\n'
	EneCG_deriv = 0
	pb = sim.utility.ProgressBar(Text = 'Processing CG frames', Steps = NFrames_CG)
	Ind = 0
	for frame in FrameRange_CG:
	   Pos_CG = Trj_CG[frame]
	   x = 0.0
	   x = getdudrc(pos = Pos_CG, nmon = poly_dict['N_mon']*poly_dict['N_poly'],
					boxl = BoxL, ldcut = LDCut,coeff = coeff, dcoeff = dcoeff, 
					nknot = nknot, spdist = spdist, spc1 = spc1, spc2 = spc2, spc3 = spc3)
	   
	   if DEBUG: print "CG->", x
	   EneCG_deriv += x
	   
	   pb.Update(Ind)
	   Ind += 1
			
	dSrel = beta * (EneAA_deriv/NFrames_AA - EneCG_deriv/NFrames_CG)
	if DEBUG: print "dSrel->",dSrel

	# logging data
	of = open(datalog, 'a')
	of.write('%g \t %g\n' % (LDCut, dSrel))
	of.close()
	

	
