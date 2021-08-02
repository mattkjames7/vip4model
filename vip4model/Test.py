import numpy as np
from .Model import Model

def Test(R=0.85,MaxDeg=4):
	'''
	This is a simple function to test the model by recreating a plot in 
	Connerney et al 2018 (figure 4, sort of).
	
	Inputs
	======
	R : float
		The radial distance to evaluate the model at.
	MaxDeg : int
		Maximum model degree to calculate.
	
	'''
	try:
		import matplotlib.pyplot as plt
		import matplotlib.colors as colors
		from mpl_toolkits.axes_grid1 import make_axes_locatable
	except:
		raise SystemError('This function requires "matplotlib" to be instaled')

	#get the coordinates to calculate the model at
	lat = np.linspace(-90,90,181)
	lon = np.linspace(0.0,360.0,361)
	latc = 0.5*(lat[1:] + lat[:-1])
	lonc = 0.5*(lon[1:] + lon[:-1])
	long,latg = np.meshgrid(lon,lat)
	longc,latgc = np.meshgrid(lonc,latc)

	longcr = longc*np.pi/180.0
	latgcr = (90.0 - latgc)*np.pi/180.0
	r = np.zeros(longcr.shape) + R
	
	#calculate the model
	Br,Bt,Bp = Model(r,latgcr,longcr,MaxDeg)
	
	#B = np.sqrt(Br**2 + Bt**2 + Bp**2)
	
	#convert to Gauss
	Bg = Br.reshape(longcr.shape)*1e-5
	
	

	#get the scale
	scale = [-30.0,30.0]
	
	#set norm
	norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
		
	maps = [1,1,0,0]
	fig = plt
	fig.figure()
	ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	ax.set_aspect(1.0)
	ax.set_xlabel('SIII East Longitude ($^\circ$)')
	ax.set_ylabel('SIII Latitude ($^\circ$)')
		
	sm = ax.pcolormesh(long,latg,Bg,cmap='RdYlBu_r',norm=norm)
	ct = ax.contour(longc,latgc,Bg,colors='grey',levels=np.linspace(-25,25,11))
	ax.clabel(ct, inline=True, fontsize=8,fmt='%2d')

	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)

	cbar = plt.colorbar(sm,cax=cax) 
	cbar.set_label('$B_r$ (Gauss) at $r$ = {:4.2f}'.format(R) + ' R$_{j}$')
	
	return ax
