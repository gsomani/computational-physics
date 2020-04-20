import numpy as np
import matplotlib.pyplot as plt

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

G=1
M=10
L=2

def f(r,t):
	d=r[0]*r[0]+r[1]*r[1]
	R=(d+0.25*L*L)**0.5
	a=-G*M/(d*R)	    
	vx=r[2]
	vy=r[3]
	ax=a*r[0]
	ay=a*r[1]
	return np.array([vx,vy,ax,ay])

r,t=rk4([1,0,0,1],10,2**-10,f)
x=r[:,0]
y=r[:,1]
plt.plot(x,y)
plt.grid(True, which="both")
plt.gca().set_aspect('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Orbit') 
plt.savefig('trajectory.png') 
plt.show()
