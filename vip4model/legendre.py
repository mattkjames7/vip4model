import numpy as np
from .binomial import binomial

def legendre(n,x):
	
	Pn = 0.0
	for k in range(0,n+1):
		Pn += binomial(n,k)*binomial(n+k,k)*((x-1)/2.0)**k
	
	return Pn
	
def legendregrad(n,x):
	
	dPndx = 0.0
	for k in range(0,n+1):
		dPndx += binomial(n,k)*binomial(n+k,k)*(k/2)*((x-1)/2.0)**(k-1)
	
	return dPndx

def oddfact(n):
	'''
	1*3*5..(2n-1)
	'''
	out = 1
	end = 2*n -1
	i = 1
	while (i <= end):
		out = out*i
		i += 2
	return out 

def delta(i,j):
	return np.int32(i == j)

def schmidtcoeff(n,m):
	
	Snm = np.sqrt(((2-delta(m,0))*np.math.factorial(n-m) )/np.math.factorial(n+m))*(oddfact(n)/np.math.factorial(n-m))
	return Snm
