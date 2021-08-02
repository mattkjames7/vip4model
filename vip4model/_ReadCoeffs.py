import numpy as np
from . import Globals
import PyFileIO as pf

def _ReadCoeffs():
	'''
	Read in the file containing the g and h Schmidt normalized 
	coefficients.
	
	Returns
	=======
	out : numpy.recarray
		Structured array containing all of the model coefficients.
		Contains the following fields:
		n : int
			n index
		m : int
			m index
		g : float
			g coefficient
		h : float
			h coefficient

	'''
	
	#the data file
	fname = Globals.ModulePath + '/__data/coeffs.dat'
	
	#read it in
	coef = pf.ReadASCIIData(fname,Header=False,dtype=[('gh','U1'),('n','int32'),('m','int32'),('coef','float32')])
	
	#create the output object
	dtype = [	('n','int32'),
				('m','int32'),
				('g','float32'),
				('h','float32')]
	n = np.sum(np.arange(4)+2)
	out = np.recarray(n,dtype=dtype)
	out.g[:] = 0.0
	out.h[:] = 0.0
	
	p = 0
	for n in range(1,5):
		for m in range(0,n+1):
			out.n[p] = n
			out.m[p] = m
			p += 1
			
	#fill the array in
	for c in coef:
		ind = np.sum(np.arange(c.n)+1)-1 + c.m 
		out[c.gh][ind] = c.coef
	
	#convert grom Gauss to nT
	out.g *= 100000.0
	out.h *= 100000.0
		
	return out
