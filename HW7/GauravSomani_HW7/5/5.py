import random
import numpy as np

def volume(m):
	if(m==1):
		return 2
	if(m==2):
		return np.pi
	return 2*np.pi*volume(m-2)/m

def volume_mc(m,n):
	total=0
	for i in range(n):
		cur = 0
		for j in range(m):
			cur += random.random()**2
		total += cur<1
	I = total/n
	error = (I*(1-I)/n)**0.5 
	return 2**m*I,2**m*error

m=10
volume_10 = volume(m)
vol10_mc,err = volume_mc(m,4**(m+1))
diff = abs(vol10_mc - volume_10)

print("Volume of %i-dimensional hypersphere = %.2f, Estimated error = %.2f, Actual error = %.2f" % (m,vol10_mc,diff,err))
