import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def tdma_reduced(a,b,c):
    N=len(b)
    for i in range(N):
        b[i] -= c[i-1]*a[i]
        c[i] /= b[i]
    return b,c

def solve_tdma_reduced(a,b,c,d):
    N=len(d)
    for i in range(N):
        d[i]=(d[i]-d[i-1]*a[i])/b[i]
    for i in range(N-2,-1,-1):
        d[i] -= c[i]*d[i+1]
    return d

hbar = constants.hbar
m = 9.109e-31
L = 1e-8
x0,sigma,k = L/2 , 1e-10, 5e10
period=4*m*L*L/(hbar*np.pi)

def f0(x):
    return np.exp(1j*k*x-(x-x0)*(x-x0)/(2*sigma*sigma))

def init(dt,N):
    psi=np.zeros(N+1,dtype=complex)
    x = np.linspace(0,L,N+1)
    psi[1:N] = f(x[1:N])
    psi[0]=psi[N]=0    
    dx = L/N
    alpha = 1j*dt*hbar/(4*m*dx*dx)
    a2,a1 = (1+2*alpha)*np.ones(N-1),-alpha*np.ones(N-1)
    a3 = 1*a1
    a1[0],a3[N-2] = 0,0    
    b1,b2 = 1-2*alpha,alpha          
    a2_,a3_=tdma_reduced(a1,a2,a3)
    return x,psi,a1,a2_,a3_,b1,b2 

def update(psi,a1,a2,a3,b1,b2):    
    N=len(psi)-1
    y=np.zeros(N-1,dtype=complex)
    for i in range(1,N):
        y[i-1] = b1*psi[i]+b2*(psi[i+1]+psi[i-1])
    psi[1:N] = solve_tdma_reduced(a1,a2,a3,y)
    return psi

dt = 2**-59
N = 1000
time = 2**-53
t=0

x,psi,a1,a2,a3,b1,b2=init(dt,N)

fig,ax = plt.subplots()
line, = ax.plot(psi)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    line.set_ydata(np.abs(psi))
    return line,time_text

def animate(i,psi):
    plt.grid(True,'both')
    psi=update(psi,a1,a2,a3,b1,b2)
    ax.set_ylabel('Re $\psi(x)$')
    ax.set_xlabel('x (m)')
    ax.set_title('Gaussian wavepacket in a box')
    line.set_ydata(np.abs(psi))
    time_text.set_text('time = %.2Es' % (i*dt))
    if(i%32==0):
        plt.savefig('wavefunction'+str(int(i/32))+'.png')
    return line,time_text

ani = animation.FuncAnimation(fig, animate, frames=2**11, fargs=(psi,),interval=100, blit=True,repeat = False)
plt.show()
