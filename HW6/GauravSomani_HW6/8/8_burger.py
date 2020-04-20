import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

a = 1
u0_deriv_min = -a*2*np.pi
tb = -1/u0_deriv_min

def f(u):
    return 0.5*u*u

def u0(x):
    return a*np.sin(2*np.pi*x)   

def u(N):
    x = np.linspace(0,1,N+1)
    return x,u0(x)

def lax(u,f,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        u_[i] = 0.5*(u[i+1] + u[i-1] - dt*N*(f(u[i+1]) - f(u[i-1])))
    return u_

def upwind_conservative(u,f,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        if u[i]>0:
            diff = - dt*(f(u[i]) - f(u[i-1]))*N
        elif u[i]<0:
            diff = - dt*(f(u[i+1]) - f(u[i]))*N
        else:
            diff = - 0.5*dt*(f(u[i+1]) - f(u[i-1]))*N
        u_[i] = u[i] + diff
    return u_

def upwind(u,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        if u[i]>0:
            diff = - dt*u[i]*(u[i] - u[i-1])*N
        elif u[i]<0:
            diff = - dt*u[i]*(u[i+1] - u[i])*N
        else:
            diff = 0
        u_[i] = u[i] + diff
    return u_

def lax_wendroff(u,f,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        u_[i] = 0.5*(u[i+1] + u[i] - dt*(f(u[i+1]) - f(u[i]))*N)
    for i in range(N):
        u[i] -= dt*(f(u_[i]) - f(u_[i-1]))*N
    return u

def lax_wendroff_viscous(u,f,l,dt):   
    N=len(u)
    qcon = l*N
    u_= np.empty(N)
    def D_coeff(u):
        D=np.zeros(N)
        for i in range(-1,N-1):
            diff=u[i+1]-u[i-1]
            if(diff<0):
                D[i]=l*l*N*abs(diff)*N*N
        return D
    dt_vis = 1/(4*qcon*N)
    M = round(dt/dt_vis)
    for j in range(M//2):
        D = D_coeff(u)
        for i in range(-1,N-1):
            u_[i] = u[i] + D[i]*dt_vis*(u[i+1] - 2*u[i] + u[i-1])
        u = np.copy(u_)
    for i in range(-1,N-1):
        u_[i] = 0.5*(u[i+1] + u[i] - dt*(f(u[i+1]) - f(u[i]))*N)
    for i in range(N):
        u[i] -= dt*(f(u_[i]) - f(u_[i-1]))*N
    for j in range(M//2):
        D = D_coeff(u)
        for i in range(-1,N-1):
            u_[i] = u[i] + D[i]*dt_vis*(u[i+1] - 2*u[i] + u[i-1])
        u = np.copy(u_)
    return u

N=2**7
x,U = u(N)
dt = 0.5/N
t=0
u_lax,u_upwind,u_upwind_conservative,u_lax_wendroff=np.copy(U),np.copy(U),np.copy(U),np.copy(U)
time = [0,0.1,0.25,0.5]

for i in range(4):
    while(t<time[i]):
        t+=dt
        u_lax[:N]=lax(u_lax[:N],f,dt)
        u_upwind[:N]=upwind(u_upwind[:N],dt)
        u_upwind_conservative[:N]=upwind_conservative(u_upwind_conservative[:N],f,dt)
        u_lax_wendroff[:N]=lax_wendroff(u_lax_wendroff[:N],f,dt)
    u_lax[N],u_upwind[N],u_upwind_conservative[N],u_lax_wendroff[N]=u_lax[0],u_upwind[0],u_upwind_conservative[0],u_lax_wendroff[0]
    plt.plot(x,u_lax,label='Lax')
    plt.plot(x,u_upwind,label='Upwind')
    plt.plot(x,u_upwind_conservative,label='Upwind (Conservative)')
    plt.plot(x,u_lax_wendroff,label='Lax-Wendroff')
    plt.ylim([-1.1,1.1])
    plt.ylabel('u')
    plt.xlabel('x')
    plt.title('{Burger\'s equation} t = '+str(time[i])+' (128 grid points)')
    plt.grid(True,'both')
    plt.gca().legend()
    plt.savefig('u128_'+str(i)+'.png')
    plt.show()

t=0
N=2**8
x,u_lax_wendroff_v = u(N)
for i in range(4):
    while(t<time[i]):
        t+=dt
        u_lax_wendroff_v[:N]=lax_wendroff_viscous(u_lax_wendroff_v[:N],f,2/N,dt)
    u_lax_wendroff_v[N]=u_lax_wendroff_v[0]
    plt.plot(x,u_lax_wendroff_v)
    plt.ylim([-1.1,1.1])
    plt.ylabel('u')
    plt.xlabel('x')
    plt.title('Burger\'s equation solution at t = '+str(time[i])+'(Lax-Wendroff with artificial viscosity )')
    plt.grid(True,'both')
    plt.savefig('artificial_viscosity'+str(i)+'.png')
    plt.show()
