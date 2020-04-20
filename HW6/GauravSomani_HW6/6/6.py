import numpy as np
import matplotlib.pyplot as plt
from dcst import dst,idst
from scipy import constants
import matplotlib.animation as animation

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

pi = np.pi
hbar = constants.hbar
m = 9.109e-31
L = 1e-8
x0,sigma,k = L/2 , 1e-10, 5e10
N = 1000
w0 = hbar*pi*pi/(2*m*L*L)

def f(x):
    return np.exp(1j*k*x-(x-x0)*(x-x0)/(2*sigma*sigma))

def init(N):
    psi=np.zeros(N,dtype=complex)
    x = np.linspace(0,L,N+1)[:N]
    psi = f(x)
    psi[0] = 0
    return np.real(psi),np.imag(psi)

ps_r,ps_i=init(N)
alpha,eta = dst(ps_r),dst(ps_i)
plt.show()

def psi_r(t,alpha,eta):
    a=np.zeros(N)
    for i in range(1,N):
        a[i]=alpha[i]*np.cos(i*i*w0*t)+eta[i]*np.sin(i*i*w0*t)
    ps = np.zeros(N+1)    
    ps[:N] = idst(a)
    return ps,np.linspace(0,L,N+1)

period=4*m*L*L/(hbar*np.pi)
time=2**-53
dt=2**-59

psi,x=psi_r(time,alpha,eta)
plt.plot(x,psi)
plt.grid(True,'both')
plt.ylabel('Re $\psi(x)$')
plt.xlabel('x (m)')
plt.title('Gaussian wavepacket in a box (t = $10^{-16}$ s)')
plt.savefig('1e16_wavepacket.png')
plt.show()

fig,ax = plt.subplots()
line, = ax.plot(psi)
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    line.set_ydata(np.real(psi))
    return line,time_text

def animate(i):
    psi,x = psi_r(i*dt,alpha,eta)
    line.set_ydata(psi)
    plt.grid(True,'both')
    ax.set_ylabel('Re $\psi(x)$')
    ax.set_xlabel('x (m)')
    ax.set_title('Gaussian wavepacket in a box')
    time_text.set_text('time = %.2Es' % (i*dt))
    if(i%32==0):
        plt.savefig('wavefunction'+str(int(i/32))+'.png')
    return line,time_text

ani = animation.FuncAnimation(fig, animate, frames=2**11,interval=1,blit=True, repeat = False)
plt.show()
