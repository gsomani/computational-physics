import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

pi=np.pi
year=365
D = 0.1
a,b,tau = 10,12,year
L = 20
T_bot=11

def T0(t):
    return a + b*np.sin(2*pi*t/tau)

def init(Ti,N):
    T=Ti*np.ones(N+1)
    T[0]=T0(0)
    T[N]=T_bot
    return T

def update(m,t,T):    
    N=len(T)-1
    diff=np.zeros(N+1)
    for i in range(1,N):
        diff[i]=m*(T[i+1]+T[i-1] - 2*T[i])
    T=T+diff
    T[0]=T0(t+h)
    return T,t+h

def solve(D,h,N,L,Ti):
    m=D*h*(N/L)**2
    x=np.linspace(0,L,N+1)
    T=init(10,N)
    t=0
    while(t<9*year):
        T,t=update(m,t,T)
    t=0
    prof = []
    for i in range(1,5):
        while(4*t < i*year):
            T,t=update(m,t,T)
        prof.append(T)
    return prof,x

h=2**-4
N=160
Ti=10
T,x=solve(D,h,N,L,Ti)
for i in range(4):
    plt.plot(x,T[i],label=str(i*3)+' months')
plt.grid(True,'both')
plt.xlabel('Depth (m)')
plt.ylabel('Temperature ($^{\circ}C$)')
plt.title('Temperature profile')
plt.gca().legend()
plt.savefig('temperature.png')
plt.show()
