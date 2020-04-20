import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import math

def forward_euler(x0,time,h,A):
    d=len(x0)
    N=int(round(1+time/h))
    x,t=np.empty([N,d]),np.empty(N)
    x[0],t[0]=x0,0
    B=np.transpose(h*A)
    for i in range(1,N):
        x[i],t[i] = x[i-1] + np.dot(x[i-1],B),t[i-1]+h
    return x,t

def backward_euler(x0,time,h,A):
    d=len(x0)
    N=int(round(1+time/h)) 
    x,t=np.empty([N,d]),np.empty(N)
    x[0],t[0]=x0,0
    I=np.identity(d)
    B=linalg.inv(I-h*A)
    B=np.transpose(B)
    for i in range(1,N):
        x[i],t[i] = np.dot(x[i-1],B),t[i-1]+h 
    return x,t

A = np.matrix('998 1998; -999 -1999') 
e = np.array([-1,-1000])
t_stab = 2/max(abs(e))
C = np.matrix('2 -1; -1 1')

def U(t):
    f= [math.exp(e[0]*t),math.exp(e[1]*t)]
    return np.dot(C,f)

print("Maximum step size allowed for stable explicit forward euler scheme = %.3f" %t_stab)

time=50*t_stab

N=1000
y=np.empty([N+1,2])
T=np.empty([N+1])

for i in range(N+1):
    T[i]=time*i/N
    y[i]=U(T[i])
   
u=[1,0]
var=['u','v']
h=[0.5*t_stab,10*t_stab]
title=['\Delta t_{stab}/2','10\Delta t_{stab}']
st=['stable','unstable']

for j in range(2):
    u_f,t=forward_euler(u,time,h[j],A)
    u_b,t=backward_euler(u,time,h[j],A)
    for i in range(2):
        plt.grid(True, which="both")
        plt.plot(t,u_f[:,i],label="explicit")
        plt.plot(t,u_b[:,i],label="implicit")
        plt.plot(T,y[:,i],label="analytic")
        plt.ylim([-0.05+min(y[:,i]),max(y[:,i])+0.05])
        plt.xlabel('t')
        plt.ylabel(var[i]+'(t)')
        plt.title(var[i]+'(t) {'+'$\Delta t = '+title[j]+'$}') 
        plt.gca().legend()
        plt.savefig(var[i]+'_'+st[j]+'.png')
        plt.show()
