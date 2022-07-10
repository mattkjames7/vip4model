import numpy as np
from ._SphHarm import _SphHarm,_SphHarmScalar

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
	_r = np.array([r]).flatten()*1.0023695021241394
	_t = np.array([theta]).flatten()
	_p = np.array([phi]).flatten()
	
	return _SphHarm(_r,_t,_p,MaxDeg)
	
def ModelScalar(r,theta,phi,MaxDeg=10):
	'''
	JRM09 Magnetic field model (see Connerney et al 2018 below). The 
	model uses right-handed System III coordinates (I think). 
	
	Inputs
	======
	r : float
		Radial distance in Rj.
	theta : float
		Colatitude in radians.
	phi : float
		East longitude in radians.
	MaxDeg : int
		Maximum degree of the model, valid 1 - 10 (default = 10).
		
	Returns
	=======
	Br : float
		Radial field
	Bt : float
		Meridional field
	Bp : float
		Azimuthal field
		
	If using the JRM09 model, please cite the following paper:
	
	Connerney, J. E. P., Kotsiaros, S., Oliversen, R. J., Espley, J. R., 
	Joergensen, J. L., Joergensen, P. S., et al. (2018). A new model of 
	Jupiter's magnetic field from Juno's first nine orbits. Geophysical 
	Research Letters, 45, 2590– 2596. https://doi.org/10.1002/2018GL077312
	
	'''
	

	return _SphHarmScalar(r*1.0023695021241394,theta,phi,MaxDeg)
