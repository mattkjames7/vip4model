import numpy as np
from ._SphHarm import _SphHarm

def Model(r,theta,phi,MaxDeg=4):
	'''
	VIP4 magnetic field model (see Connerney et al 1998 below). The 
	model uses right-handed System III coordinates (I think). 
	
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
		Meridional field
	Bp : float
		Azimuthal field
	
	Please cite:
	
	Connerney, J. E. P., Acuña, M. H., Ness, N. F., and Satoh, T. (1998),
	New models of Jupiter's magnetic field constrained by the Io flux 
	tube footprint, J. Geophys. Res., 103( A6), 11929– 11939, 
	doi:10.1029/97JA03726.
		
	'''
	
	#make all the inputs arrays
	_r = np.array([r]).flatten()
	_t = np.array([theta]).flatten()
	_p = np.array([phi]).flatten()
	
	return _SphHarm(_r,_t,_p,MaxDeg)
