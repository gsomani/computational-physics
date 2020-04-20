import scipy as sp
from scipy import constants
import math
import matplotlib.pyplot as plt

c=constants.c
h=constants.h
hbar=constants.hbar
pi=constants.pi
k_b=constants.k
sigma_unicode="\u03EC"

def E(x):
    if(x==0):
        return 0
    return x**3/(math.exp(x)-1)

def g(z):
    if(z==1):
        return 0
    x=z*3.5/(1-z)
    return 3.5*E(x)/((z-1)*(z-1))

N=200
x=sp.linspace(0,1,N+1)
total=0

for i in range(1,N-1,2):
	total+=4*g(x[i])+2*g(x[i+1])
total+=g(x[0])+4*g(x[N-1])+g(x[N])
total/=3*N

sigma=(k_b**4)*total/(c*c*h*h*hbar)

print(sigma_unicode+"=%.3E (in SI units)" %sigma)
