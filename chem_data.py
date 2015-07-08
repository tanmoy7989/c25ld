#!/usr/bin/env python
import numpy as np

def getModelData():
	""" Returns dictionaries containing data about the polymer-water system
	    and about the different cutoffs used in the system
	"""
	poly_dict = {"N_mon": 25, "N_poly": 1,  "N_water" : 1700,
             	     "blength"      :   1.5300,     "Mass_mon"  :  16.0427,
              	     "Kbond"        :   400,   	     "r0"       :   1.53,
             	     "Kangle"       :  55.17204,    
             	     "theta0"    : 111.00*(np.pi/180.0),
             	     "LJEpsilon"    :  0.13986,     "LJSigma"   :    3.73}

	cut_dict = {"SPCut": 7.0, "LJCut" : 10.0, "WCACut" : 4.187,
				"LDCut": 6.5, "LDLowerCut" : 4.635}
				
	TempSet = 298.0
	ret = {'poly_dict': poly_dict, 'cut_dict': cut_dict, 'TempSet': TempSet}
	return ret
        
