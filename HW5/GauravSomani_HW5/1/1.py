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
    N=int(round(1+time/h)) 
    x,t=np.empty(N),np.empty(N)
    x[0],t[0]=x0,0
    for i in range(1,N):
        x[i],t[i]=rk4_step(x[i-1],t[i-1],h,f)
    return x,t

def vin(t):
    x=int(2*t)
    if(x%2==0):
        return 1
    return -1

def f(vout,t):
    return (vin(t)-vout)/rc

p=[0.01,0.1,1]

plt.grid(True, which="both")

for i in range(3):
    rc=p[i]
    x,t=rk4(0,10,rc/(2**6),f)
    plt.plot(t,x,label='RC = %.2f s' % rc)

plt.xlabel('t (in seconds)')
plt.ylabel('$V_{out}$') 
plt.title('Low pass filter')
plt.xticks(np.linspace(0,10,11))
plt.gca().legend()
plt.savefig('filter.png')
plt.show()
            
