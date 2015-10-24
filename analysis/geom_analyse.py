#!/usr/bin/env python

import os, sys
import numpy as np
import pickle
import sim
import pickleTraj #kept in $HOME/.scripts

## GLOBALS
Measures = ['Rg', 'R_EE', 'SASA_atom', 'kappa', 'pmf_Rg_R_EE']
NMeasures = len(Measures)

def compute(trj, trajtype, measure, nmon = 25, npoly = 1, nwater = 1700):
    
    doMinImage = True
    
    cgtype = trajtype[1]
    if cgtype == 'AA':  
        poly_indices = range(3*nwater, 3*nwater+nmon)
        start = 0; stop = len(trj); freq = 100
    else:   
        poly_indices = range(0,nmon)
        start = 0; stop = len(trj); freq = 10
    
    Nframes = int((stop - start)/freq)
    ind = 0
    if measure == 'SASA_atom':
        output = np.zeros([Nframes, nmon], np.float64)
    else:
        output = np.zeros(Nframes)
    
    pb = sim.utility.ProgressBar(Text = 'Processing frames....', Steps = Nframes)
    for frame in xrange(start, stop, freq):
        Pos = trj[frame][poly_indices]
        if doMinImage:
             for i in range(1,nmon):
                Pos[i] = sim.geom.Reimage(Pos[i], Pos[i-1], trj.FrameData['BoxL'])
       
        if measure == 'Rg':
            output[ind] = np.sqrt(np.sum((Pos-Pos.mean(axis=0))**2)/nmon) 
        elif measure == 'R_EE':
            output[ind] = np.sqrt(np.sum((Pos[-1]-Pos[0])**2))
        elif measure == 'kappa':
            l = np.zeros(3)
            diff = Pos - Pos.mean(axis = 0)
            for ax in range(3): l[ax] = np.sqrt(np.sum(diff[:,ax]**2, axis = 0))
            output[ind] = ((1.5) * (l[0]**4 + l[1]**4 + l[2]**4) / (l[0]**2 + l[1]**2 + l[2]**2)**2) - 0.5
        elif measure == 'SASA_atom':
            mon_rad = 2.0934 ; w_rad = 1.4
            radii = (mon_rad + w_rad) * np.ones(nmon)
            output[ind, :] = sim.geom.SphereSurfaceAreas(Pos, radii)
        
        ind += 1    
        pb.Update(ind) 
    
    if measure == 'SASA_atom': output = np.ndarray.flatten(output)    
    return output   
    
   
def makeHist1D(z, nbins = 50, normalize = True):
    Nframes = len(z)
    bin_min = 0.98*z.min() ; bin_max = 1.02*z.max()
    delta = (bin_max - bin_min)/float(nbins)
    bin_centers = np.zeros(nbins); bin_vals = np.zeros([nbins], np.float64)
    for i in range(nbins):
        bin_centers[i] = bin_min + (i+0.5)*delta
    
    for i in range(Nframes):
        assignment = int((z[i] - bin_min)/delta)
        bin_vals[assignment] += 1.0
    
    bin_vals /= Nframes
    if normalize: 
        norm = np.trapz(bin_vals, bin_centers, dx = delta)
        bin_vals /= norm
  
    return (bin_centers, bin_vals)
    
    
def makeHist2D(z0, z1, nbins = (125, 125), normalize = True):
    if not len(z0) == len(z1):
        raise TypeError('The 2 measures must have same # of computed samples')
    
    Nframes = len(z0)        
    bin_min = (0.98*z0.min(), 0.98*z1.min())
    bin_max = (1.02*z0.max(), 1.02*z1.max())
    delta = ((bin_max[0]-bin_min[0])/nbins[0], (bin_max[1]-bin_min[1])/nbins[1])
    bin_centers_0 = np.zeros(nbins[0]); bin_centers_1 = np.zeros(nbins[1])
    bin_vals = np.zeros([nbins[0], nbins[1]], np.float64)

    for i in range(nbins[0]): bin_centers_0[i] = z0.min() + (i+0.5)*delta[0]
    for i in range(nbins[1]): bin_centers_1[i] = z1.min() + (i+0.5)*delta[1]
    
    for i in range(Nframes):
        assignment_0 = int((z0[i]-bin_min[0])/delta[0])
        assignment_1 = int((z1[i]-bin_min[1])/delta[1])
        bin_vals[assignment_0, assignment_1] += 1.0
        
    #bin_vals /= Nframes
    if normalize:
        Row_norm = np.zeros(nbins[1])
        for i in range(nbins[1]):
            Row_norm[i] = np.trapz(bin_vals[i,:], bin_centers_1, dx = delta[1])
        
        #norm = np.trapz(Row_norm, bin_centers_0, dx = delta[0])
        norm = np.sum(bin_vals) * delta[0] * delta[1]
        bin_vals /= norm
    
    return ((bin_centers_0, bin_centers_1), bin_vals)
    
    
def writeToFile(prefix, data):
    filename = prefix + '.pickle'
    pickle.dump(data, open(filename, 'w'))


def readFromFile(prefix):
    filename = prefix + '.pickle'
    return pickle.load(open(filename, 'r'))

        
def getFilePrefix(basedir, data_name, data_type = 'measure', trajtype = None):
    dirname = os.path.join(basedir, data_type)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    
    filePrefix = os.path.join(dirname, '%s_%s_%s' % (data_name, trajtype[0], trajtype[1]))
    return filePrefix

                
def isComputed(basedir, measure, trajtype):
    prefix = getFilePrefix(basedir = basedir, data_name = measure, data_type = 'measure', trajtype = trajtype)
    if os.path.isfile(prefix+'.pickle'): return True
    else: return False
    
    
def isBinned(basedir, hist, trajtype):
    prefix = getFilePrefix(basedir = basedir, data_name = hist, data_type = 'hist', trajtype = trajtype)
    if os.path.isfile(prefix+'.pickle'): return True
    else: return False
    
       
def main(TrajFile, trajtype, N_mon = 25, N_poly = 1, N_water = 1700, savedir = os.getcwd()):

    Trj = pickleTraj(os.path.abspath(TrajFile))
    
    # Computing
    for measure in Measures:
        if measure == 'pmf_Rg_R_EE' or isComputed(savedir, measure, trajtype):    continue  
        print 'Computing %s...' % measure
        data = compute(trj = Trj, trajtype = trajtype, measure = measure, nmon = N_mon, npoly = N_poly, nwater = N_water)
        prefix = getFilePrefix(basedir = savedir, data_name = measure, data_type = 'measure', trajtype = trajtype)
        writeToFile(prefix, data)
    
    # Making histograms    
    for measure in Measures:
        if isBinned(savedir, measure, trajtype):  continue
        print 'Binning %s...' % measure
        hist_prefix = getFilePrefix(basedir = savedir, data_name = measure, data_type = 'hist', trajtype = trajtype)
        
        if measure != 'pmf_Rg_R_EE':
            compute_prefix = getFilePrefix(basedir = savedir, data_name = measure, data_type = 'measure', trajtype = trajtype)
            data = readFromFile(compute_prefix)
            hist = makeHist1D(data)
            writeToFile(hist_prefix, hist)
        
        if measure == 'pmf_Rg_R_EE':
            Rg_prefix = getFilePrefix(basedir = savedir, data_name = 'Rg', data_type = 'measure', trajtype = trajtype)
            R_EE_prefix = getFilePrefix(savedir, data_name = 'R_EE', data_type = 'measure', trajtype = trajtype)
            
            Rg = readFromFile(Rg_prefix); R_EE = readFromFile(R_EE_prefix)
            hist = makeHist2D(Rg, R_EE)
            writeToFile(hist_prefix, hist)



### User inputs
TrajFile = sys.argv[1]

trajtype_string = sys.argv[2]
fftype = trajtype_string.split(',')[0][1:]
cgtype = trajtype_string.split(',')[1][:-1]
trajtype = (fftype, cgtype)

N_mon = int(sys.argv[3]); N_poly = int(sys.argv[4]) ; N_water = int(sys.argv[5])
savedir = sys.argv[6]

main(TrajFile, trajtype, N_mon, N_poly, N_water, savedir)                        
