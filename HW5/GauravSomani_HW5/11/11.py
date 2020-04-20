import numpy as np
from scipy import constants
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

V0=50
a=1e-11
hbar=constants.hbar
e=constants.e
m=constants.m_e
gamma=2*m*e*V0*a*a/(hbar*hbar)

n=3
eps=2**-20
h=2**-7
xr=[-10*a,10*a]

def rk4_step(g,x,h,E,n,f):
    k0=f(g,x,n,E)
    k1=f(g+0.5*h*k0,x+0.5*h,n,E)
    k2=f(g+0.5*h*k1,x+0.5*h,n,E)
    k3=f(g+h*k2,x+h,n,E)
    gf = g + h*(k0+2*k1+2*k2+k3)/6
    return gf,x+h

def rk4(g0,xr,h,f,n,E):
    N=int(round(1+(xr[1]-xr[0])/h))
    g,x=np.empty([N,2]),np.empty(N)
    g[0],x[0]=g0,xr[0]
    for i in range(1,N):
        g[i],x[i]=rk4_step(g[i-1],x[i-1],h,E,n,f)
    return g,x

def f(r,x,n,E):
    f_psi=r[1]
    f_derivative=gamma*(x**(2*n) - E)*r[0]
    return np.array([f_psi,f_derivative])

def secant_root(x,f,eps):
    while(abs(x[1]-x[0])>=eps/2):
        y=[f(x[0]),f(x[1])]
        m=(y[1]-y[0])/(x[1]-x[0])
        x = x[1],x[1]-y[1]/m
    return x[1]

def normalising_factor(psi,xr,norm):
    total=0
    length=xr[1]-xr[0]
    N=len(psi)-1
    g=np.multiply(psi,psi)
    for i in range(1,N-1,2):
        total+=4*g[i]+2*g[i+1]
    total+=g[0]+4*g[N-1]+g[N]
    total*=length/(3*N)
    return (a*total/norm)**0.5

def phi(E,xr,n):
    phi0=[0,2**-52]
    g,x=rk4(phi0,xr,h,f,n,E)
    return g,x

def solution(Er,xr,n,eps,m):
    xr = xr[0]/a,xr[1]/a
    def boundary_value(E):
        g,x=phi(E,xr,n)
        y=g[:,0]
        z=g[:,1]
        N=len(y)
        if(m%2==0):
            return z[N-1]
        return y[N-1]
    
    Er = Er[0]/V0,Er[1]/V0
    E=secant_root(Er,boundary_value,eps)
    y,x_half=phi(E,xr,n)
    g=y[:,0]
    N=len(g)
    g=g/normalising_factor(g,[xr[0],0],0.5)
    psi=np.empty(2*N-1)
    psi[0:N]=g
    x=np.empty(2*N-1)
    x[0:N]=x_half   
    if(m%2==0):
        psi[N-1:2*N-1]=g[::-1]
    else:
        psi[N-1:2*N-1]=-g[::-1]
    x[N-1:2*N-1]=-x_half[::-1]
    return psi,a*x,E*V0

E=np.empty([n])
osc=["Harmonic","Anharmonic"]

xr=[-10*a,0]

for j in range(2):
    k=j+1
    Ec=gamma**(-k/(k+1))*V0
    E0=Ec
    psi,x,E[0]=solution([E0,E0+Ec],xr,k,eps,0)
    N=len(x)
    plt.plot(x,psi,label='n=0')
    E0=3*E[0]
    print("For "+osc[j]+" oscillator,")
    print('Ground state energy = %.3f eV'%(E[0]))
    for i in range(1,n):
        psi,x,E[i]=solution([E0,E0+Ec],xr,k,eps,i)
        E0=2*E[i]-E[i-1]
        plt.plot(x,psi,label='n='+str(i))
        diff=E[i]-E[i-1]
        print('E = %.3f eV for n = %i ; E_%i - E_%i = %.3f eV'%(E[i],i,i,i-1,diff))
    print('\n')
    plt.gca().legend()
    plt.grid(True,"both")
    if(j==1):
        xlimit=[-5*a,5*a]
    else:
        xlimit=[-10*a,10*a]
    plt.xlim(xlimit)
    plt.xticks(np.linspace(xlimit[0],xlimit[1],11))
    plt.xlabel('x (m)')
    plt.ylabel('$\psi(x)$')
    plt.title("V(x)=$V_0 (x/a)^"+str(2*k)+"$")
    plt.savefig(str(osc[j])+'_oscillator.png')
    plt.show()
