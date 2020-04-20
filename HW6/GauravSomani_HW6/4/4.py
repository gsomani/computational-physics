import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

N=100
T=10
g=9.8
       
def relaxation(f,h,w):
    err=0
    for i in range(1,N):
        diff = w*((f[i-1]+f[i+1])/2 - f[i] + g*(h*h/2))
        f[i]+=diff
        err=max(abs(diff),err)
    return f,err

eps=2**-20   
    
def solve(eps,N):
    x = np.zeros([N+1])
    err=eps
    w=2/(1+np.sin(np.pi/N))
    h=T/N
    while(err>eps/2):
        x,err=relaxation(x,h,w)
    return x,np.linspace(0,T,N+1)

def exact_solution(N):
    x = np.zeros([N+1])
    H=g*T*T/8
    t=np.linspace(0,T,N+1)
    x[N//2:N]=H-0.5*g*(t[:N//2])**2    
    x[N//2:0:-1]=x[N//2:N]
    return x,t

y,t=exact_solution(N)
x,t=solve(eps,N)
plt.plot(t,x)
plt.grid(True,'both')
plt.xlabel('Time (s)')
plt.ylabel('Height (m)')
plt.title('Height of ball as a function of time')
plt.savefig('height.png')
plt.show()
