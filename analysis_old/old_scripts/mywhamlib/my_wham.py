import numpy as np
import sys
sys.dont_write_bytecode = True

# Fortran routine written by MSS
import whamlib

kB = 1.3806e-23 * 6.023e23 / 1000.0  #Boltzmann constant in kJ/mol
Temp = 298.0
beta = 1./(kB*Temp)

def main(beta, eng, umbrella, hist, step):
	
	# Unwrapping dicts passed to the function
	nstages = step['nstages']
	nframes = step['nframes']
	
	E_pe = eng['pe']
	E_u  = eng['u']
	Rg = umbrella['Rg']
	
	bin_assignment = hist['bin_assignment']
	bin_center = hist['bin_center']
	nbins = hist['nbins']
	count = hist['count']

	Rg0 = umbrella['Rg0']
	FConst = umbrella['FConst']
	
	pmf = whamlib.free_energy(ekn = E_pe,betak = beta*np.ones(nstages),
							 nbin = nbins,niterbin = 50,niterall = 100)
							 
	return (bin_center, pmf)
	
	
	
	
	
