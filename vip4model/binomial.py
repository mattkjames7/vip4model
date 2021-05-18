import numpy as np


def _binomial(n,k) -> int:
	'''
	Calculate the binomial coefficient.
	'''
	
	if (k < 0) or (k > n):
		return 0
	elif (k == 0) or (k == n):
		return 1
	
	#symmetry
	k = np.min([k,n-k])
	c = 1
	for i in range(0,k):
		c *= ((n - i)/(i + 1))
	return c


def binomial(n,k):
	
	scalar = True
	if hasattr(n,'__iter__'):
		nn = n
		kk = k
		scalar = False
	else:
		nn = np.array([n])
		kk = np.array([k])
		
	out = np.zeros(np.size(nn),dtype='float64')
	for i in range(0,np.size(nn)):
		out[i] = _binomial(nn[i],kk[i])
		
	if scalar:
		return out[0]
	else:
		return out
