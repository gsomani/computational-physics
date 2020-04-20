import scipy as sp
import scipy.constants as constants
import math
import matplotlib.pyplot as plt

m=9.1094e-31
e=constants.e
hbar=constants.hbar
w=1e-9
alpha=2*m*e*w*w/(hbar*hbar)
V=20
phi=alpha*V

def f(theta):
    return (2*theta*theta-phi)*math.sin(theta)-2*theta*math.cos(theta)*(phi-theta*theta)**0.5

def false_position(bracket,eps):
    a,b=bracket
    fa=f(a)
    fb=f(b)
    c=(a*fb-b*fa)/(fb-fa)
    if(min(c-a,b-c) < 0.5*eps): 
        return c
    fc=f(c)
    if((fa<0 and fc>0) or (fa>0 and fc<0)):
        return false_position([a,c],eps)
    return false_position([c,b],eps)


pi=constants.pi
b=[[0.5*pi,pi]]
theta=[]
E=[]
eps=1e-3
e=0.5*eps*(alpha/V)**0.5

for i in range(7):
    theta.append(false_position(b[i],e))
    b.append([0.5*pi+theta[i],(i+2)*pi])
    E.append(round((theta[i]*theta[i]/alpha),4))

print("Energy levels (in eV) = " + str(E))
