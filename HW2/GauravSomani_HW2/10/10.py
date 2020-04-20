import scipy as sp
from scipy import constants
import math
import cmath
import matplotlib.pyplot as plt

pi=constants.pi

def f(z):
    return cmath.exp(2*z)   

def g(m):
    N=10000
    total=0
    j=complex(0,1)
    for i in range(N):         
        z=cmath.exp(j*2*pi*i/N)
        total+=f(z)*cmath.exp(-j*2*pi*i*m/N)
    return math.factorial(m)*total/N

print("1st derivative = %.2f" %abs(g(1)) )        
print("2nd derivative = %.2f" %abs(g(2)) )        

for i in range(3,21):
    print("%ith derivative = %.2f" %(i,abs(g(i))) )
