import numpy as np
from . import Globals
from ._ReadCoeffs import _ReadCoeffs

def _CoeffGrids():
	
	#check if we have loaded them yet
	if not Globals.g is None and not Globals.h is None:
		return Globals.g,Globals.h
	
	
	#read the file with the coefficients in first
	if Globals.coeffs is None:
		Globals.coeffs = _ReadCoeffs()
	coeffs = Globals.coeffs
	
	#create the grids (n,n) shape
	MaxDeg = np.max(coeffs.n)
	g = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64')
	h = np.zeros((MaxDeg+1,MaxDeg+1),dtype='float64')
		
	#fill it in
	g[coeffs.n,coeffs.m] = coeffs.g
	h[coeffs.n,coeffs.m] = coeffs.h

	#add to globals
	Globals.g = g
	Globals.h = h
	
	return g,h
	
	
