import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
import time

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

G = 6.6738e-20
M = 1.9891e30
GM = G*M
distance=4e9
velocity=0.5
day=86400
year=365*day
delta=1/year

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
        t[i] = theta[i] + e*np.sin(theta[i])
    x[1+N//2:],y[1+N//2:] = x[N//2-1:0:-1],-y[N//2-1:0:-1]
    t[1+N//2:] = 2*np.pi - t[N//2-1:0:-1]
    t *= time_per_radian
    x=x+distance
    return a,b,e,T,x,y,t

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
    pos = x[:,0:2]
    vel = x[:,2:4]
    return pos,t

def adaptive_rk4(x0,time,h,delta,f,mult):
    x,t=x0,0
    X=[x0[0]]
    Y=[x0[1]]
    T=[t]
    h0=h
    i=1
    mean_h=0
    step=[]
    eps=1/4
    while(time>t):
        xh,th = rk4_step(x,t,h,f)
        x1,tf = rk4_step(xh,th,h,f)
        x2,tf = rk4_step(x,t,2*h,f)
        diff = x1 - x2
        d=la.norm(diff[:2])
        m=2*(h*delta/d)**(1/4)
        if(m < 1):
            h *= m 
            continue
        t=tf 
        T.append(th)
        T.append(t)
        mean_h+=h
        step.append(h)
        x = x1 + diff/15
        X.append(xh[0])	
        Y.append(xh[1])
        X.append(x[0])	
        Y.append(x[1])
        diff_t=time-t
        if(diff_t<eps):
            break	
        h*=min(mult,m)
        i+=2
        if(diff_t<2*h):
            h=diff_t/2
    mean_h /= (i-1)//2
    min_h,max_h=min(step),max(step)
    return X,Y,T,min_h,mean_h,max_h

def f(r,t):
	d=-GM/la.norm(r[:2])**3
	vx=r[2]
	vy=r[3]
	ax=d*r[0]
	ay=d*r[1]
	return np.array([vx,vy,ax,ay])

N=2**17
a,b,e,T,xe,ye,te=trajectory(distance,velocity)

h=T/N
start_time = time.time()
r,t=rk4([distance,0,0,velocity],2*T,h,f)
time_fixed=time.time()-start_time
x,y = r[:,0],r[:,1]
drift_fixed = 0.5*year*la.norm([x[N-1]-distance,y[N-1]])/T
print("Fixed time step :")
print("h = %.2E year" %(h/year))
print("Computational time (fixed time-step) = %.1f s" % time_fixed)
print("Average drift in trajectory = %.1f km/year" % drift_fixed)
plt.plot(x,y)
plt.gca().set_aspect('equal')
plt.grid(True, which="both")
plt.xlabel('x (in km)')
plt.ylabel('y (in km)')
plt.title('Orbit (Fixed time step)')
plt.savefig('trajectory_fixed.png') 
plt.show()   

start_time = time.time()
x,y,t,min_h,h,max_h=adaptive_rk4([distance,0,0,velocity],2*T,1.5*year,delta,f,1.1)
time_adaptive=time.time()-start_time
drift_adaptive = 0.5*year*la.norm([x[len(t)-1]-distance,y[len(t)-1]])/T
print("\nAdaptive time step :")
print("Minimum h = %.1E year, Mean h = %.1E year, Maximum h = %.1f year" % (min_h/year,h/year,max_h/year))
print("Computational time (adaptive time-step) = %.1f s" % time_adaptive)
print("Average drift in trajectory = %.2f km/year" % drift_adaptive)
ratio_drift=drift_fixed/drift_adaptive
ratio_time=time_fixed/time_adaptive
print("Average drift in trajectory for adaptive time-step method is smaller compared to fixed time-step by a factor of %.1E " % ratio_drift)
print("Computational time taken by adaptive method is smaller by a factor of %i as compared to fixed time-step method" % ratio_time)
print("\nSo, adaptive method is faster and more accurate compared to fixed time-step method.")
plt.plot(x,y,'orange')
plt.plot(x[:-1+len(x)//2:2],y[:-1+len(y)//2:2],'.',markersize=2,mec='blue')
plt.gca().set_aspect('equal')
plt.grid(True, which="both")
plt.xlabel('x (in km)')
plt.ylabel('y (in km)')
plt.title('Orbit (Adaptive time step)')
plt.savefig('trajectory_adaptive.png') 
plt.show()
