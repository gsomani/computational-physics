import numpy as np
import matplotlib.pyplot as plt

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

def f(r,t):
    fx=a*r[0]-b*r[0]*r[1]
    fy=-d*r[1]+g*r[0]*r[1]
    return np.array([fx,fy])
   
a=1
b=g=0.5
d=2

r,t=rk4([2,2],30,2**-10,f)
x=r[:,0]
y=r[:,1]
plt.plot(t,x,label='Rabbits')
plt.plot(t,y,label='Foxes')
plt.grid(True, which="both")
plt.xlabel('t')
plt.ylabel('Population (in thousands)')
plt.title('Population dynamics of rabbits and foxes') 
plt.gca().legend()
plt.savefig('prey.png')
plt.show()
            
