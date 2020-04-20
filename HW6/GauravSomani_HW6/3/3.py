import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

C,d,sigma,L,v=1,0.1,0.3,1,100

def f(x):
    return C*x*(L-x)*np.exp(-(x-d)*(x-d)/(2*sigma*sigma))/(L*L)

def init(N):
    e=np.zeros(N+1)
    u=np.zeros(N+1)
    x=np.linspace(0,L,N+1)
    u[1:N] = f(x[1:N])
    return e,u,x

def update(e,u,h,N):
    e_=np.copy(e)
    u_=np.copy(u)
    m=v*v*h*N*N/(L*L)
    e_[1:N] += h*u[1:N]
    for i in range(1,N):
        u_[i] += m*(e[i+1]+e[i-1]-2*e[i])
    return e_,u_

N=1000
h=2**-20
e,u,x=init(N)
stable=True
i=1
while(stable):
    e,u=update(e,u,h,N)
    if(i%200==0):
        stable = np.min(e)>=0
        plt.plot(x,1e3*e,label='t = '+str(round(i*1000*h,2))+' ms')
    i+=1
plt.xlabel('x (in m)')
plt.ylabel('$\\xi(x)$ {in mm}')
plt.grid(True,'both')
plt.title('Transverse displacement {$\\xi(x)$}')
plt.gca().legend()
plt.savefig('displacement.png')
plt.show()
