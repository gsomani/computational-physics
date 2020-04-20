import random
import numpy as np

def random_gen():
	return random.random()**2

def g(x):
	return 1/(np.exp(x)+1)

def integral(g,n):
	I=0
	varf = 0
	for i in range(n):
		x = random_gen()
		I += g(x)
	return 2*I/n

n = 2**20

area = integral(g,n)
print("Integral = %.2f" % area)
