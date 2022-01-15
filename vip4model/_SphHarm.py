import numpy as np
from ._CoeffGrids import _CoeffGrids
from ._Legendre import _Legendre,_LegendreScalar
from ._Schmidt import _Schmidt

def _SphHarm(r,theta,phi,MaxDeg=4):
	'''
	This function calculates the VIP4 model field using spherical
	harmonics. 
	
	Inputs
	======
	r : float
		Radial coordinate in Rj (69,911 km).
	theta : float
		Colatitude in radians (RH SIII).
	phi : float
		Azimuth in radians (RH SIII).
	MaxDeg : int
		Maximum degree of the model to use - can be 1 to 4 (default = 4)
		for VIP4, where larger degrees take longer to process, but are
		better.
	
	Returns
	=======
	Br : float
		Magnetic field strength in radial direction (nT).
	Bt : float
		Magnetic field strength in latitudinal direction (nT).
	Bp : float
		Magnetic field strength in azimuthal direction (nT).
	
	'''

	#calculate cosmphi and sinmphi
	nr = np.size(r)
	cosmphi = np.zeros((MaxDeg+1,nr),dtype='float64') + 1.0
	sinmphi = np.zeros((MaxDeg+1,nr),dtype='float64')
	for m in range(1,MaxDeg+1):
		mphi = m*phi
		cosmphi[m] = np.cos(mphi)
		sinmphi[m] = np.sin(mphi)
	
	#get the grids of g and h parameters
	g,h = _CoeffGrids()
	
	#output arrays
	Br = np.zeros(r.size,dtype='float64')
	Bt = np.zeros(r.size,dtype='float64')
	Bp = np.zeros(r.size,dtype='float64')
	

	#start calculating p costheta and its derivative
	Pnm,dPnm = _Legendre(theta,MaxDeg)

	#get the Schmidt coefficients
	#not sure if these are needed - I wonder if the g and h coefficients
	#are already normalised...
	Snm = _Schmidt(MaxDeg)
	
	#Normalize the polynomials
	SP = (Pnm.T*Snm.T).T
	SdP = (dPnm.T*Snm.T).T

	#temporary arrays for the inner sums
	sumr = np.zeros((r.size,),dtype='float64')
	sumt = np.zeros((r.size,),dtype='float64')
	sump = np.zeros((r.size,),dtype='float64')
	
	#now sum everything up
	r1 = 1/r
	sintheta1 = 1.0/np.sin(theta)
	badst = np.where(np.isfinite(sintheta1) == False)[0]
	if badst.size > 0:
		sintheta1[badst] = 0.0
	C = (r1)**2
	for n in range(1,MaxDeg+1):
		#define the constant (a/r)**(n+2)
		C = C*r1

		#do the innder sum for this value of n first
		#I suspect this could be vectorised, but would probably be 
		#pretty unreadable if I did
		sumr[:] = 0
		sumt[:] = 0
		sump[:] = 0
		for m in range(0,n+1):
			sumr += SP[n,m]*(g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m])
			sumt += SdP[n,m]*(g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m])
			sump += m*SP[n,m]*(h[n,m]*cosmphi[m] - g[n,m]*sinmphi[m])

		#add it to the output arrays
		Br += C*(n+1)*sumr
		Bt += -C*sumt
		Bp += -C*sump
	Bp = sintheta1*Bp

	return Br,Bt,Bp



def _SphHarmScalar(r,theta,phi,MaxDeg=4):
	'''
	This function calculates the VIP4 model field using spherical
	harmonics. 
	
	Inputs
	======
	r : float
		Radial coordinate in Rj (69,911 km - this is wrong).
	theta : float
		Colatitude in radians (RH SIII).
	phi : float
		Azimuth in radians (RH SIII).
	MaxDeg : int
		Maximum degree of the model to use - can be 1 to 4 (default=4)
		for VIP4, where larger degrees take longer to process, but are
		better.
	
	Returns
	=======
	Br : float
		Magnetic field strength in radial direction (nT).
	Bt : float
		Magnetic field strength in latitudinal direction (nT).
	Bp : float
		Magnetic field strength in azimuthal direction (nT).
	
	'''
	

	#calculate cosmphi and sinmphi
	cosmphi = np.zeros((MaxDeg+1),dtype='float64') + 1.0
	sinmphi = np.zeros((MaxDeg+1),dtype='float64')
	for m in range(1,MaxDeg+1):
		mphi = m*phi
		cosmphi[m] = np.cos(mphi)
		sinmphi[m] = np.sin(mphi)
	
	#get the grids of g and h parameters
	g,h = _CoeffGrids()
	
	#output arrays
	Br = 0.0
	Bt = 0.0
	Bp = 0.0
	

	#start calculating p costheta and its derivative
	Pnm,dPnm = _LegendreScalar(theta,MaxDeg)

	#get the Schmidt coefficients
	#not sure if these are needed - I wonder if the g and h coefficients
	#are already normalised...
	Snm = _Schmidt(MaxDeg)
	
	#Normalize the polynomials
	SP = Pnm*Snm
	SdP = dPnm*Snm

	#temporary arrays for the inner sums
	sumr = 0.0
	sumt = 0.0
	sump = 0.0
	
	#now sum everything up
	r1 = 1/r
	sintheta1 = 1.0/np.sin(theta)
	if np.isfinite(sintheta1) == False:
		sintheta1 = 0.0
	C = r1*r1
	for n in range(1,MaxDeg+1):
		#define the constant (a/r)**(n+2)
		C = C*r1

		#do the innder sum for this value of n first
		#I suspect this could be vectorised, but would probably be 
		#pretty unreadable if I did
		sumr = 0.0
		sumt = 0.0
		sump = 0.0
		for m in range(0,n+1):
			gchs = g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m]
			sumr += SP[n,m]*gchs
			sumt += SdP[n,m]*gchs
			sump += m*SP[n,m]*(h[n,m]*cosmphi[m] - g[n,m]*sinmphi[m])

		#add it to the output arrays
		Br += C*(n+1)*sumr
		Bt += -C*sumt
		Bp += -C*sump
	Bp = sintheta1*Bp

	return Br,Bt,Bp

def _SphHarmScalarCart(x,y,z,MaxDeg=4):
	'''
	This function calculates the VIP4 model field using spherical
	harmonics in  Cartesian coordinates. 
	
	Inputs
	======
	x : float
		Radial coordinate in Rj
	y : float
		Colatitude  in Rj
	z : float
		Azimuth  in Rj
	MaxDeg : int
		Maximum degree of the model to use - can be 1 to 4 (default=4)
		for VIP4, where larger degrees take longer to process, but are
		better.
	
	Returns
	=======
	Bx : float
		Magnetic field strength in radial direction (nT).
	By : float
		Magnetic field strength in latitudinal direction (nT).
	Bz : float
		Magnetic field strength in azimuthal direction (nT).
	
	'''
	#calculate r,theta,phi and rho
	rho = np.sqrt(x**2 + y**2)
	r = np.sqrt(x**2 + y**2 + z**2)
	theta = np.arccos(z/r)
	phi = (np.arctan2(y,x) + (2*np.pi)) % (2*np.pi)
		

	#calculate cosmphi and sinmphi
	cosmphi = np.zeros((MaxDeg+1),dtype='float64') + 1.0
	sinmphi = np.zeros((MaxDeg+1),dtype='float64')
	for m in range(1,MaxDeg+1):
		mphi = m*phi
		cosmphi[m] = np.cos(mphi)
		sinmphi[m] = np.sin(mphi)
	
	#get the grids of g and h parameters
	g,h = _CoeffGrids()
	
	#output arrays
	Bx = 0.0
	By = 0.0
	Bz = 0.0
	

	#start calculating p costheta and its derivative
	Pnm,dPnm = _LegendreScalar(theta,MaxDeg)

	#get the Schmidt coefficients
	#not sure if these are needed - I wonder if the g and h coefficients
	#are already normalised...
	Snm = _Schmidt(MaxDeg)
	
	#Normalize the polynomials
	SP = Pnm*Snm
	SdP = dPnm*Snm


	#some constants for equations 175-177
	#175
	x_o_r = x/r
	xz_o_pr2 = (x*z)/(rho*r*r)
	y_o_p = y/(rho*rho)
	#176
	y_o_r = y/r
	yz_o_pr2 = (y*z)/(rho*r*r)
	x_o_p = x/(rho*rho)
	#177
	z_o_r = z/r
	p_o_r2 = -rho/(r*r)
	
	if x == 0:
		x_o_p = 0.0
		xz_o_pr2 = 0.0
	if y == 0:
		y_o_p = 0.0
		yz_o_pr2 = 0.0
	
	print(x_o_r,xz_o_pr2,y_o_p)
	print(y_o_r,yz_o_pr2,x_o_p)
	print(z_o_r,p_o_r2)
	
	#now sum everything up
	r1 = 1/r
	C = r1*r1
	for n in range(1,MaxDeg+1):
		#define the constants (a/r)**(n+1) and (a/r)**(n+2)
		Cprev = C
		C = C*r1
		
		#n + 1
		np1 = n + 1

		#this method requires 2 inner sums for each component because
		#I am not intelligent enough to combine them!
		sumxyz0 = 0.0
		sumx1 = 0.0
		sumy1 = 0.0
		sumz1 = 0.0
		for m in range(0,n+1):
			gchs = g[n,m]*cosmphi[m] + h[n,m]*sinmphi[m]
			gshc = g[n,m]*sinmphi[m] - h[n,m]*cosmphi[m]
			sumxyz0 += SP[n,m]*gchs

			sumx1 += xz_o_pr2*SdP[n,m]*gchs + m*y_o_p*SP[n,m]*gshc
			sumy1 += yz_o_pr2*SdP[n,m]*gchs - m*x_o_p*SP[n,m]*gshc
			sumz1 += p_o_r2*SdP[n,m]*gchs

		#add it to the output arrays
		Bx += np1*x_o_r*C*sumxyz0 - Cprev*sumx1
		By += np1*y_o_r*C*sumxyz0 - Cprev*sumy1
		Bz += np1*z_o_r*C*sumxyz0 - Cprev*sumz1


	return Bx,By,Bz
