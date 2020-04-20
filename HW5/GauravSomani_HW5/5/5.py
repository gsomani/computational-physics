import numpy as np
import matplotlib.pyplot as plt
from math import pi,cos,sin

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

def rk4_step(x,t,h,f,m):
	k0=f(x,t,m)
	k1=f(x+0.5*h*k0,t+0.5*h,m)
	k2=f(x+0.5*h*k1,t+0.5*h,m)
	k3=f(x+h*k2,t+h,m)
	xf = x+h*(k0+2*k1+2*k2+k3)/6
	return xf,t+h

def rk4(x0,h,f,m):
    x=x0
    t=0
    time=2*x0[3]/g
    N=round(1+time/h) 
    X,T=np.empty([N,4]),np.empty(N)
    X[0],T[0]=x0,0
    i=1
    while(i<N):
        x,t=rk4_step(X[i-1],T[i-1],h,f,m)
        if(x[2]<0):
            h=h/2
            continue
        X[i],T[i]=x,t
        i+=1
    return X[:i],T[:i]

R=0.08
C=0.47
rho=1.22
g=9.8
drag=pi*R*R*rho*C*0.5
u=100

def f(r,t,m):
	a=drag/m
	fx=r[1]
	fvx=-a*r[1]*(r[1]*r[1]+r[3]*r[3])**0.5
	fy=r[3]
	fvy=-g-a*r[3]*(r[1]*r[1]+r[3]*r[3])**0.5
	return np.array([fx,fvx,fy,fvy])

def trajectory(th,m):
    ux = u*cos(th)
    uy = u*sin(th)
    r,t=rk4([0,ux,0,uy],2**-10,f,m)
    x=r[:,0]
    y=r[:,2]
    return x,y,x[len(x)-1]

x,y,xr=trajectory(pi/6,1)
plt.plot(x,y)
plt.grid(True, which="both")
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.gca().set_aspect('equal')
plt.title('Projectile motion')
plt.savefig('projectile.png') 
plt.show()

for m in range(1,11):
    x,y,xr=trajectory(pi/6,m)
    plt.plot(x,y,label="m = "+str(m)+" kg")
plt.grid(True, which="both")
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.gca().set_aspect('equal')
plt.title('Comparision of trajectories of canonballs of different mass') 
plt.gca().legend()
plt.savefig('comparision.png')
plt.show()
