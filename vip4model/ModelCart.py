import numpy as np
from .Model import Model

def ModelCart(x,y,z,MaxDeg=4):
	
	
	#convert to spherical polar coords
	r = np.sqrt(x**2 + y**2 + z**2)
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
	
