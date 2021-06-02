import numpy as np


def _Legendre(theta,MaxDeg=4):
	
	
	Pnm = np.zeros((MaxDeg+1,MaxDeg+1,theta.size),dtype='float64') + np.nan
	dPnm = np.zeros((MaxDeg+1,MaxDeg+1,theta.size),dtype='float64') + np.nan
	sintheta = np.sin(theta)
	costheta = np.cos(theta)
	Pnm[0,0,:] = 1.0
	Pnm[1,0,:] = costheta
	Pnm[1,1,:] = sintheta
	dPnm[0,0,:] = 0.0
	dPnm[1,0,:] = -sintheta
	dPnm[1,1,:] = costheta
	for n in range(2,MaxDeg+1):
		n21 = 2*n - 1
		for m in range(0,n+1):
			if n == m:
				Pnm[n,m] = n21*sintheta*Pnm[n-1,m-1]
				dPnm[n,m] = n21*(costheta*Pnm[n-1,m-1] + sintheta*dPnm[n-1,m-1])
			elif m == (n-1):
				Pnm[n,m] = n21*sintheta*Pnm[n-1,m-1]
				dPnm[n,m] = n21*(costheta*Pnm[n-1,m-1] + sintheta*dPnm[n-1,m-1])
			else:
				Pnm[n,m] = (1.0/(n-m))*(costheta*n21*Pnm[n-1,m] - (n + m -1)*Pnm[n-2,m])
				dPnm[n,m] = (1.0/(n-m))*(n21*(costheta*dPnm[n-1,m] - sintheta*Pnm[n-1,m]) - (n + m -1)*dPnm[n-2,m])
				

	return Pnm,dPnm
