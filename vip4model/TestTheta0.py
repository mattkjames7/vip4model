import numpy as np
import matplotlib.pyplot as plt
from ._SphHarm import _SphHarmScalar,_SphHarmScalarCart

def TestTheta0(fig=None,maps=[1,1,0,0],n=101,trange=[-0.1,0.1],R=5,MaxDeg=4):
	'''
	Test both models to see what happens near theta = 0
	
	
	'''
	
	#use two arrays of coords:
	#first set - moving along x (phi = 0 or pi)
	theta = np.linspace(trange[0],trange[1],n)
	
	theta0 = np.abs(theta)
	phi0 = np.zeros(n,dtype='float64')
	phi0[theta < 0] = np.pi
	r = np.zeros(n,dtype='float64') + R
	
	#second set along y (phi = +/- pi/2)
	theta1 = theta0
	phi1 = phi0 - np.pi/2.0
	
	#calculate their Cartesian equivalents
	cost0 = np.cos(theta0)
	sint0 = np.sin(theta0)
	cosp0 = np.cos(phi0)
	sinp0 = np.sin(phi0)
	cost1 = np.cos(theta1)
	sint1 = np.sin(theta1)
	cosp1 = np.cos(phi1)
	sinp1 = np.sin(phi1)
	x0 = r*sint0*cosp0
	y0 = r*sint0*sinp0
	z0 = r*cost0
	x1 = r*sint1*cosp1
	y1 = r*sint1*sinp1
	z1 = r*cost1
	

	
	#call the spherical harmonic functions
	Br0s = np.zeros(n,dtype='float64')
	Bt0s = np.zeros(n,dtype='float64')
	Bp0s = np.zeros(n,dtype='float64')
	Br1s = np.zeros(n,dtype='float64')
	Bt1s = np.zeros(n,dtype='float64')
	Bp1s = np.zeros(n,dtype='float64')
	Bx0c = np.zeros(n,dtype='float64')
	By0c = np.zeros(n,dtype='float64')
	Bz0c = np.zeros(n,dtype='float64')
	Bx1c = np.zeros(n,dtype='float64')
	By1c = np.zeros(n,dtype='float64')
	Bz1c = np.zeros(n,dtype='float64')
	for i in range(0,n):
		Br0s[i],Bt0s[i],Bp0s[i] = _SphHarmScalar(r[i],theta0[i],phi0[i],MaxDeg)
		Br1s[i],Bt1s[i],Bp1s[i] = _SphHarmScalar(r[i],theta1[i],phi1[i],MaxDeg)
		Bx0c[i],By0c[i],Bz0c[i] = _SphHarmScalarCart(x0[i],y0[i],z0[i],MaxDeg)
		Bx1c[i],By1c[i],Bz1c[i] = _SphHarmScalarCart(x1[i],y1[i],z1[i],MaxDeg)
	
	#convert coordinate systems
	
	#sph to cart
	Bx0s = Br0s*sint0*cosp0 + Bt0s*cost0*cosp0 - Bp0s*sinp0
	By0s = Br0s*sint0*sinp0 + Bt0s*cost0*sinp0 + Bp0s*cosp0
	Bz0s = Br0s*cost0 - Bt0s*sint0	
	Bx1s = Br1s*sint1*cosp1 + Bt1s*cost1*cosp1 - Bp1s*sinp1
	By1s = Br1s*sint1*sinp1 + Bt1s*cost1*sinp1 + Bp1s*cosp1
	Bz1s = Br1s*cost1 - Bt1s*sint1	
	
	#cart to sph
	Br0c = Bx0c*sint0*cosp0 + By0c*sint0*sinp0 + Bz0c*cost0
	Bt0c = Bx0c*cost0*cosp0 + By0c*cost0*sinp0 - Bz0c*sint0
	Bp0c =-Bx0c*sinp0       + By0c*cosp0
	Br1c = Bx1c*sint1*cosp1 + By1c*sint1*sinp1 + Bz1c*cost1
	Bt1c = Bx1c*cost1*cosp1 + By1c*cost1*sinp1 - Bz1c*sint1
	Bp1c =-Bx1c*sinp1       + By1c*cosp1
	
	
	plt.figure(figsize=(8,11))
	
	ax0 = plt.subplot2grid((3,2),(0,0))
	ax1 = plt.subplot2grid((3,2),(1,0))
	ax2 = plt.subplot2grid((3,2),(2,0))
	ax3 = plt.subplot2grid((3,2),(0,1))
	ax4 = plt.subplot2grid((3,2),(1,1))
	ax5 = plt.subplot2grid((3,2),(2,1))

	ax0.plot(theta,Br0c,label='Cartesian (along $x$)')
	ax0.plot(theta,Br1c,label='Cartesian (along $y$)')
	ax0.plot(theta,Br0s,label='Spherical (along $x$)',linestyle='--')
	ax0.plot(theta,Br1s,label='Spherical (along $y$)',linestyle='--')

	ax1.plot(theta,Bt0c,label='Cartesian (along $x$)')
	ax1.plot(theta,Bt1c,label='Cartesian (along $y$)')
	ax1.plot(theta,Bt0s,label='Spherical (along $x$)',linestyle='--')
	ax1.plot(theta,Bt1s,label='Spherical (along $y$)',linestyle='--')

	ax2.plot(theta,Bp0c,label='Cartesian (along $x$)')
	ax2.plot(theta,Bp1c,label='Cartesian (along $y$)')
	ax2.plot(theta,Bp0s,label='Spherical (along $x$)',linestyle='--')
	ax2.plot(theta,Bp1s,label='Spherical (along $y$)',linestyle='--')

	ax3.plot(theta,Bx0c,label='Cartesian (along $x$)')
	ax3.plot(theta,Bx1c,label='Cartesian (along $y$)')
	ax3.plot(theta,Bx0s,label='Spherical (along $x$)',linestyle='--')
	ax3.plot(theta,Bx1s,label='Spherical (along $y$)',linestyle='--')

	ax4.plot(theta,By0c,label='Cartesian (along $x$)')
	ax4.plot(theta,By1c,label='Cartesian (along $y$)')
	ax4.plot(theta,By0s,label='Spherical (along $x$)',linestyle='--')
	ax4.plot(theta,By1s,label='Spherical (along $y$)',linestyle='--')

	ax5.plot(theta,Bz0c,label='Cartesian (along $x$)')
	ax5.plot(theta,Bz1c,label='Cartesian (along $y$)')
	ax5.plot(theta,Bz0s,label='Spherical (along $x$)',linestyle='--')
	ax5.plot(theta,Bz1s,label='Spherical (along $y$)',linestyle='--')

	ax0.set_ylabel(r'$B_r$ (nT)')
	ax1.set_ylabel(r'$B_{\theta}$ (nT)')
	ax2.set_ylabel(r'$B_{\phi}$ (nT)')
	ax3.set_ylabel(r'$B_x$ (nT)')
	ax4.set_ylabel(r'$B_y$ (nT)')
	ax5.set_ylabel(r'$B_z$ (nT)')

	ax0.set_xlabel(r'$\theta$ ($^{\circ}$)')
	ax1.set_xlabel(r'$\theta$ ($^{\circ}$)')
	ax2.set_xlabel(r'$\theta$ ($^{\circ}$)')
	ax3.set_xlabel(r'$\theta$ ($^{\circ}$)')
	ax4.set_xlabel(r'$\theta$ ($^{\circ}$)')
	ax5.set_xlabel(r'$\theta$ ($^{\circ}$)')

	ax0.legend()
	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()
	ax5.legend()
	
	ax0.set_title('Sph Code $\mathbf{B}=$('+'{:5.0f}'.format(Bx0s[n//2])+'$\hat{x}$,'+'{:5.0f}'.format(By0s[n//2])+'$\hat{y}$,'+'{:5.0f}'.format(Bz0s[n//2])+'$\hat{z}$)')
	ax3.set_title('Cart Code $\mathbf{B}=$('+'{:5.0f}'.format(Bx0c[n//2])+'$\hat{x}$,'+'{:5.0f}'.format(By0c[n//2])+'$\hat{y}$,'+'{:5.0f}'.format(Bz0c[n//2])+'$\hat{z}$)')
	
	
	plt.tight_layout()
