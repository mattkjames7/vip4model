import numpy as np
from ._CoeffGrids import _CoeffGrids
from ._Legendre import _Legendre
from ._Schmidt import _Schmidt

def _SphHarm(r,theta,phi,MaxDeg=4):


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
