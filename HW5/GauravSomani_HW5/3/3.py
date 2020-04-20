import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 24})
plt.rcParams["figure.figsize"] = (20,20)

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

def f(R,t):
    fx=s*(R[1]-R[0])
    fy=r*R[0]-R[1]-R[0]*R[2]
    fz=R[0]*R[1]-b*R[2]
    return np.array([fx,fy,fz])
   
s=10
r=28
b=8/3

R,t=rk4([0,1,0],50,2**-10,f)
x=R[:,0]
y=R[:,1]
z=R[:,2]

plt.plot(t,y)
plt.grid(True, which="both")
plt.title('y(t)')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.savefig('y_trajectory.png')
plt.show()

plt.plot(x,z)
plt.xlabel('x')
plt.ylabel('z')
plt.title('Attractor')
plt.gca().set_aspect('equal')
plt.grid(True, which="both")
plt.savefig('attractor.png')
plt.show()
