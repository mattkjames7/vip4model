import numpy as np
from .Model import Model,ModelScalar
from ._SphHarm import _SphHarm,_SphHarmScalarCart


def ModelCart(x,y,z,MaxDeg=4):
	'''
	VIP4 Magnetic field model (see Connerney et al 1998 below). The 
	model uses right-handed System III coordinates (I think).  
	
	Inputs
	======
	x : float
		x-coordinate in Rj, R-H System III.
	y : float
		y-coordinate in Rj, R-H System III.
	z : float
		z-coordinate in Rj, R-H System III.
		
	Returns
	=======
	Bx : float
		x component of magnetic field, nT.
	By : float
		y component of magnetic field, nT.
	Bz : float
		z component of magnetic field, nT.
		
	Please cite:
		
	Connerney, J. E. P., Acuña, M. H., Ness, N. F., and Satoh, T. (1998), 
	New models of Jupiter's magnetic field constrained by the Io flux 
	tube footprint, J. Geophys. Res., 103( A6), 11929– 11939, 
	doi:10.1029/97JA03726.
	
	'''		
	
	#convert to spherical polar coords
	r = np.sqrt(x**2 + y**2 + z**2)*1.0023695021241394
	theta = np.arccos(z/r)
	phi = (np.arctan2(y,x) + (2*np.pi)) % (2*np.pi)
	
	#call the model
	Br,Bt,Bp = Model(r,theta,phi,MaxDeg)

	#convert to Cartesian (hopefully correctly...)
	cost = np.cos(theta)
	sint = np.sin(theta)
	cosp = np.cos(phi)
	sinp = np.sin(phi)
	Bx = Br*sint*cosp + Bt*cost*cosp - Bp*sinp
	By = Br*sint*sinp + Bt*cost*sinp + Bp*cosp
	Bz = Br*cost - Bt*sint
	
	return Bx,By,Bz
	

def ModelCartScalar(x,y,z,MaxDeg=4):
	'''
	VIP4 Magnetic field model. The 
	model uses right-handed System III coordinates (I think).  
	
	Inputs
	======
	x : float
		x-coordinate in Rj, R-H System III.
	y : float
		y-coordinate in Rj, R-H System III.
	z : float
		z-coordinate in Rj, R-H System III.
		
	Returns
	=======
	Bx : float
		x component of magnetic field, nT.
	By : float
		y component of magnetic field, nT.
	Bz : float
		z component of magnetic field, nT.
		
	If using the VIP4 model, please cite the following paper:
	
	Connerney, J. E. P., Acuña, M. H., Ness, N. F., and Satoh, T. (1998), 
	New models of Jupiter's magnetic field constrained by the Io flux 
	tube footprint, J. Geophys. Res., 103( A6), 11929– 11939, 
	doi:10.1029/97JA03726.12
	
	'''	
	
	#convert to spherical polar coords
	r = np.sqrt(x**2 + y**2 + z**2)*1.0023695021241394
	theta = np.arccos(z/r)
	phi = (np.arctan2(y,x) + (2*np.pi)) % (2*np.pi)
	
	#call the model
	Br,Bt,Bp = ModelScalar(r,theta,phi,MaxDeg)

	#convert to Cartesian (hopefully correctly...)
	cost = np.cos(theta)
	sint = np.sin(theta)
	cosp = np.cos(phi)
	sinp = np.sin(phi)
	Bx = Br*sint*cosp + Bt*cost*cosp - Bp*sinp
	By = Br*sint*sinp + Bt*cost*sinp + Bp*cosp
	Bz = Br*cost - Bt*sint
	
	return Bx,By,Bz
	

def ModelTest(x,y,z,MaxDeg=4):
	
	
	return _SphHarmScalarCart(x,y,z,MaxDeg)
