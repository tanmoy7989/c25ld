#!/usr/bin/env python

import numpy as np
import sys
import matplotlib.pyplot as plt

import alg


def testbennett():
    kT = 0.6
    beta = 1.
    np.random.seed(12345)
    x01 = 0.5; x02 = 0.5
    alpha1 = 7.0; alpha2 = 5.0
    
    Nparams = 5
    alpha1_list = [1, 1.2, 100., 20., 200., 0.01]
    alpha2_list = [1.02, 2., 200., 21., 201., 1.02]
    x01_list = [5., 6., 7., 8., 1., 10., 200.]
    x02_list = [5., 6., 7., 8., 1., 10., 200.]

    fe_exact = np.zeros(Nparams)
    fe_bennett = np.zeros(Nparams)

    for n in range(0,Nparams):
        alpha1 = alpha1_list[n]
        alpha2 = alpha2_list[n]
        x01 = x01_list[n]
        x02 = x02_list[n]
        
        sigma1 = 1./np.sqrt(2*alpha1*beta) ; sigma2 = 1./np.sqrt(2*alpha2*beta)
        
        Nframes = 1e4
        x1 = np.random.normal(x01, sigma1, Nframes)
        x2 = np.random.normal(x02, sigma2, Nframes)
    
        BE11 =  beta * alpha1*(x1 - x01)**2.
        BE22 =  beta * alpha2*(x2 - x02)**2.
        BE12 =  beta * alpha2*(x1 - x02)**2.        
        BE21 =  beta * alpha1*(x2 - x01)**2.
    
        fe1, err = alg.BennettFE(BE11 = BE11, BE22 = BE22, BE12 = BE12, BE21 = BE21, Verbose = False)
        fe2 = - 0.5 * np.log(alpha1/alpha2)
    
        fe_exact[n] = fe2
        fe_bennett[n] = fe1
        
    print'\n\n-----------------------------\n'
    print 'Case \t Exact \t Bennett\n'
    for n in range(Nparams):
        print '%d \t %g \t %g\n' % (n, fe_exact[n], fe_bennett[n])
    print'\n\n-----------------------------\n'

    
def testphi():
    LDCut = 6.5
    Delta = 1.2
    r = np.linspace(LDCut-Delta-0.5, LDCut+2.0, 100)
    coeff, dcoeff = alg.calcCoeff(LDCut = LDCut, Delta = Delta)
    phi = np.zeros(100)
    
    for n in range(100):
        if r[n] <= LDCut - Delta: phi[n] = 1.0
        elif r[n] >= LDCut: phi[n] = 0.0
        else:   phi[n] = coeff[0] + coeff[1]*r[n]**2. + coeff[2]*r[n]**4. + coeff[3]*r[n]**6.
    plt.figure
    plt.plot(r, phi, linewidth = 3)
    plt.xlabel(r'$r$'); plt.ylabel(r'$\phi(r)$')
    

def testdphidrc():
    LDCuts = [4., 5., 6.,]
    Delta = 1.2
    r = np.linspace(0,8., 100)
    
    for LDCut in LDCuts:
        coeff, dcoeff = alg.calcCoeff(LDCut = LDCut, Delta = Delta)
        dphidrc = dcoeff[0] + dcoeff[1]*r**2. + dcoeff[2]*r**4. + dcoeff[3]*r**6.
        plt.figure
        plt.plot(r, dphidrc, label = 'LDCUT = %g' % LDCut, linewidth = 3)

    plt.legend(); plt.xlabel(r'$r$'); plt.ylabel(r'$\frac{d\phi}{dr_c} (r)$')

        
def testSplineCoeff():
    import parse_potential as pp
    
    prefix = '/home/cask0/home/tsanyal/c25ld/data/cg_ff/methane/methane_wca_SPLD'
    logfile, histfile, sumfile = pp.parseFileNames(prefix)
    paramstring = pp.parseParamString(sumfile)
    h = pp.parseHist(histfile)['LD']
    
    rho_t = h[0]; rho_m = h[1] 
       
    d = alg.makeSysData(isPolymer = False, hasLJ = True)
    Sys = alg.makeSys(LDCut = 6.5, SysData = d, Prefix = 'testsys', paramstring = paramstring, fftype = 'wca')
    (nknot, spdist, spc0, spc1, spc2, spc3) = alg.getLDSpline(Sys = Sys)
    
    rho = rho_t
    dfdrho = np.zeros(len(rho))
    frho = np.zeros(len(rho))
    
    for n in range(len(rho)):
        x = (rho[n] - min(rho)) * spdist[1]
        i = max(min(int(x), nknot-2), 0)
        t = max(min(x - float(i), 1.0), 0.0)
        frho[n] = spc0[i] + t * (spc1[i] + t*(spc2[i]+t*spc3[i]))
        dfdrho[n] = spdist[1] * (spc1[i] + t * (2.*spc2[i]+t*3.*spc3[i]))
    
    fig = plt.figure
    ax1 = plt.subplot(211) ; ax2 = plt.subplot(212)
    ax1.plot(rho, frho, linewidth = 3, label = r'$f(\rho)$'); ax1.legend()
    ax2.plot(rho, dfdrho, linewidth = 3, label = r'$\frac{df}{d\rho}$'); ax2.legend()



if __name__ == '__main__':
    testSplineCoeff()
    #testdphidrc()
    #testphi()
    plt.show()
    #testbennett()
