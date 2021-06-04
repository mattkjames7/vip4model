import numpy as np
from ._SphHarm import _SphHarm

def Model(r,theta,phi,MaxDeg=4):
	'''
	
	Inputs
	======
	r : float
		Radial distance in Rj.
	theta : float
		Colatitude in radians.
	phi : float
		Azimuth in radians
		
	Returns
	=======
	Br : float
		Radial field
	Bt : float
		Meridional? field
	Bp : float
		Azimuthal field
	'''
	
	#make all the intputs arrays
	_r = np.array([r]).flatten()
	_t = np.array([theta]).flatten()
	_p = np.array([phi]).flatten()
	
	return _SphHarm(_r,_t,_p,MaxDeg)
