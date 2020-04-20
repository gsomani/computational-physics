import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def leapfrog_step(x,y,t,h,f):
	xf = x + h*f(y,t+0.5*h) 
	t=t+h
	yf = y + h*f(xf,t)
	return xf,yf,t

def leapfrog(x0,time,h,f):
    d=len(x0)
    N=int(round(1+time/h)) 
    x,t=np.empty([N,d]),np.empty(N)
    x[0],t[0]=x0,0
    y=x0+0.5*h*f(x0,0)
    for i in range(1,N):
        x[i],y,t[i]=leapfrog_step(x[i-1],y,t[i-1],h,f)
    return x,t

def f(r,t):
    fx=r[1]
    fy=r[1]*r[1]-r[0]-5
    return np.array([fx,fy])
   
r,t=leapfrog([1,0],50,1e-3,f)
x=r[:,0]
plt.plot(t,x)
plt.grid(True, which="both")
plt.xticks(np.linspace(0,50,11))
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('x(t) {Leapfrog method}')
#plt.savefig('trajectory.png') 
plt.show()
