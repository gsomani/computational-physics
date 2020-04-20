import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def rk4_step(x,t,h,f):
	k0=f(x,t)
	k1=f(x+0.5*h*k0,t+0.5*h)
	k2=f(x+0.5*h*k1,t+0.5*h)
	k3=f(x+h*k2,t+h)
	xf = x+h*(k0+2*k1+2*k2+k3)/6
	return xf,t+h

def rk4(x0,time,h,f):
    d=len(x0)
    N=int(round(1+time/h)) 
    x,t=np.empty([N,d]),np.empty(N)
    x[0],t[0]=x0,0
    for i in range(1,N):
    	x[i],t[i]=rk4_step(x[i-1],t[i-1],h,f)
    return x,t

m=1
k=6
w=2
N=5

def f(r,t):
    F=np.empty(2*N)
    Fd=math.cos(w*t)
    for i in range(1,N-1):
        F[i]=r[N+i]
        F[N+i]=(k/m)*(r[i+1]+r[i-1]-2*r[i])
    F[N]=(k/m)*(r[1]-r[0]) + Fd/m
    F[2*N-1]=(k/m)*(r[N-2]-r[N-1])
    F[0]=r[N]
    F[N-1]=r[2*N-1]    
    return F    

r_init=np.zeros(2*N)
r,t=rk4(r_init,20,2**-10,f)
for i in range(N):
    x=r[:,i]
    plt.plot(t,x,label='$\\xi_{'+str(i+1)+'}$')

plt.gca().legend()
plt.grid(True, which="both")
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Displacement of masses')
plt.savefig('displacement.png') 
plt.show()            
