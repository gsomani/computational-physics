import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

G = 6.6738e-11
M = 1.9891e30
GM=G*M
m = 5.9722e24

def trajectory(distance,velocity):
    alpha=0.5*(1+1/(2*GM/(distance*velocity*velocity)-1))
    a=distance*alpha
    e=abs(1/alpha - 1)
    b=a*(1-e*e)**0.5
    time_per_radian = a*a*(1-e*e)**0.5/(distance*velocity)
    T=2*np.pi*time_per_radian
    N=2**10
    x=np.empty(N)
    y=np.empty(N)
    t=np.empty(N)
    theta=np.linspace(0,np.pi,1+N//2)
    for i in range(1+N//2):
        x[i],y[i]=a*(np.cos(theta[i])-1),b*np.sin(theta[i])
        t[i] = theta[i] - e*np.sin(theta[i])
    x[1+N//2:],y[1+N//2:] = x[N//2-1:0:-1],-y[N//2-1:0:-1]
    t[1+N//2:] = 2*np.pi - t[N//2-1:0:-1]
    t *= time_per_radian
    x=x+distance
    return a,b,e,T,x,y,t

def verlet_step(x,v,t,h,f):
	xf = x + h*v
	pe=-GM*m/la.norm(xf)	 
	t=t+h
	vt = v + 0.5*h*f(xf,t)
	ke = 0.5*m*la.norm(vt)**2	 	
	vf = v + h*f(xf,t)
	return xf,vf,ke,pe,t

def verlet(x0,v0,time,h,f):
    N=int(round(1+time/h)) 
    x=np.empty([N,2])
    ke=np.empty(N)
    pe=np.empty(N)
    t=np.empty(N)
    x[0]=x0
    t[0]=0
    pe[0]=-GM*m/la.norm(x0)
    ke[0]=0.5*m*la.norm(v0)**2	
    v=v0+0.5*h*f(x0,0)
    for i in range(1,N):
    	x[i],v,ke[i],pe[i],t[i] = verlet_step(x[i-1],v,t[i-1],h,f)
    te=ke+pe        
    return x,ke,pe,te,t

def f(r,t):
	d=-GM/la.norm(r)**3
	ax=d*r[0]
	ay=d*r[1]
	return np.array([ax,ay])

distance=1.471e11
velocity=3.0287e4

x0=[distance,0]
v0=[0,velocity]

a,b,e,T,xe,ye,te=trajectory(distance,velocity)

hour=3600
year=365*24*hour
   
r,ke,pe,te,t=verlet(x0,v0,3*year,hour,f)
x=r[:,0]
y=r[:,1]
plt.plot(x,y)
plt.gca().set_aspect('equal')
plt.grid(True, which="both")
plt.xlabel('x (in m)')
plt.ylabel('y (in m)')
plt.title('Orbit')
plt.savefig('trajectory.png') 
plt.show()

t=t/hour
plt.plot(t,ke,label="Kinetic Energy")
plt.plot(t,pe,label="Potential Energy")
plt.plot(t,te,label="Total Energy")
plt.grid(True, which="both")
plt.xlabel('Time (in hours)')
plt.ylabel('Energy (in J)')
plt.gca().legend()
plt.title('Energy variation')
plt.savefig('energy.png') 
plt.show()

plt.plot(t,te)
plt.grid(True, which="both")
plt.xlabel('Time (in hours)')
plt.ylabel('Energy (in J)')
plt.title('Total energy variation')
plt.savefig('total_energy.png') 
plt.show()
