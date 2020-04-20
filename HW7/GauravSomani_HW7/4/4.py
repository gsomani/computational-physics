import random
import numpy as np

def hit_miss(f,xr,yr,n):
	total=0
	diff = [xr[1]-xr[0],yr[1]-yr[0]]
	for i in range(n):
		x,y = xr[0]+diff[0]*random.random(),yr[0]+diff[1]*random.random()
		total += y<f(x)
	A = diff[0]*diff[1]
	I = total*A/n
	error = (I*(A-I)/n)**0.5
	return I,error

def mean_value(f,xr,n):
	I=0
	diff = xr[1]-xr[0]
	varf = 0
	for i in range(n):
		x = xr[0]+diff*random.random()
		I += f(x)
		varf += f(x)*f(x)
	I /= n
	varf = (varf/n - I*I)
	I *= diff
	error = diff*(varf/n)**0.5
	return I,error

def f(x):
	return np.sin(1/(x*(2-x)))**2

n = 2**14

area_1,error_1 = hit_miss(f,[0,2],[0,1],n)
area_2,error_2 = mean_value(f,[0,2],n)

print("Area (hit-miss) = %.3f, Error = %.3f" %(area_1,error_1))
print("Area (mean value) = %.3f, Error = %.3f" %(area_2,error_2))
