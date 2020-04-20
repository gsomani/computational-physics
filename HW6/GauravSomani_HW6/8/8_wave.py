import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

v=1

def f(N,f0,t):
    x = np.linspace(0,1,N+1)
    x_ = x-v*t
    y = f0(x_)
    return x,y

def f0(x):
    x -= np.floor(x)
    return np.exp(-200*(x-0.3)*(x-0.3)) + (x[:]>=0.6) + (x[:]<=0.8) - 1

def lax(u,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        u_[i] = 0.5*(u[i+1] + u[i-1] - v*N*dt*(u[i+1] - u[i-1]))
    return u_

def upwind(u,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        u_[i] = u[i] + v*N*dt*(u[i-1] - u[i])
    return u_

def lax_wendroff(u,dt):
    N=len(u)
    u_= np.empty(N)
    for i in range(-1,N-1):
        u_[i] = 0.5*(u[i+1] + u[i] - v*dt*(u[i+1] - u[i])*N)
    for i in range(N):
        u[i] -= v*dt*(u_[i] - u_[i-1])*N
    return u

def lax_wendroff_v(u,l,dt):
	N=len(u)
	u_= np.empty(N)
	qcon = l*N
	def D_coeff(u):
		D=np.zeros(N)
		for i in range(-1,N-1):
			diff=u[i+1]-u[i-1]
			D[i]=l*l*N*abs(diff)*N*N
		return D
	dt_vis = 1/(4*qcon*N)
	M = int(round(dt/dt_vis))
	for j in range(M//2):
		D = D_coeff(u)
		for i in range(-1,N-1):
			u_[i] = u[i] + D[i]*dt_vis*(u[i+1] - 2*u[i] + u[i-1])
		u = np.copy(u_)
	for i in range(-1,N-1):
		u_[i] = 0.5*(u[i+1] + u[i] - v*dt*(u[i+1] - u[i])*N)
	for i in range(N):
		u[i] -= v*dt*(u_[i] - u_[i-1])*N
	for j in range(M//2):
		D = D_coeff(u)
		for i in range(-1,N-1):
			u_[i] = u[i] + D[i]*dt_vis*(u[i+1] - 2*u[i] + u[i-1])
		u = np.copy(u_)
	return u

def l1_error(f_numeric,t):
    N=len(f_numeric)    
    x,F = f(N,f0,t)
    diff = sum(abs(f_numeric-F[:N]))
    return diff/N

N=2**7
x,F = f(N,f0,0)
dt = 0.5/(v*N)
t=0

f_lax,f_upwind,f_lax_wendroff,f_lax_wendroff_v=np.copy(F),np.copy(F),np.copy(F),np.copy(F)
for i in range(4):
	while(t<i):
		t+=dt
		f_lax[:N]=lax(f_lax[:N],dt)
		f_upwind[:N]=upwind(f_upwind[:N],dt)
		f_lax_wendroff[:N]=lax_wendroff(f_lax_wendroff[:N],dt)
		f_lax_wendroff_v[:N]=lax_wendroff_v(f_lax_wendroff_v[:N],1/N,dt)
	f_lax[N],f_upwind[N],f_lax_wendroff[N]=f_lax[0],f_upwind[0],f_lax_wendroff[0]    
	plt.plot(x,F,label='Analytical')
	plt.plot(x,f_lax,label='Lax')
	plt.plot(x,f_upwind,label='Upwind')
	plt.plot(x,f_lax_wendroff,label='Lax-Wendroff')
	plt.plot(x,f_lax_wendroff_v,label='Lax-Wendroff (with artificial viscosity)')
	plt.ylabel('f')
	plt.xlabel('x')
	plt.title('{Advection equation} t = '+str(i)+' (128 grid points)')
	plt.grid(True,'both')
	plt.gca().legend()
	plt.savefig('f128_'+str(i)+'.png')
	plt.show()

N=2**5
l1=np.empty([4,4,5])
n=2**(5+np.arange(5))
dx = 1/n
slope=np.empty([4,4])

def plot_error(l1,j): 
    method = ['Lax','Upwind','Lax-Wendroff','Lax-Wendroff (with artificial viscosity)']
    for i in range(4):
        plt.loglog(dx,l1[i],label=method[i],basex=2,basey=2)
        slope[j][i]=0.5*np.log2(l1[i][2]/l1[i][4])
        print("Slope of L1-error plot at t = %.1f is %.2f for %s" %(time[j],slope[j][i],method[i]))
    plt.ylabel('L1 Error')
    plt.xlabel('$\Delta x$')
    plt.title('{Advection equation} L1 error {t = '+ str(time[j]) + '}')
    plt.grid(True,'both')
    plt.gca().legend()
    plt.savefig('l1_error_advection_'+'t_'+str(j)+'.png')
    plt.show()

time=[0.5,1,2,3]
g = []
for i in range(5):
	N=n[i]
	dt = 0.5/(v*N)
	t=0
	x,F = f(N,f0,0)
	f_lax,f_upwind,f_lax_wendroff,f_lax_wendroff_v=np.copy(F),np.copy(F),np.copy(F),np.copy(F)
	while(t<time[len(time)-1]):
		t+=dt
		f_lax[:N]=lax(f_lax[:N],dt)
		f_upwind[:N]=upwind(f_upwind[:N],dt)
		f_lax_wendroff[:N]=lax_wendroff(f_lax_wendroff[:N],dt)
		f_lax_wendroff_v[:N]=lax_wendroff_v(f_lax_wendroff_v[:N],1/N,dt)
		for j in range(len(time)):
			if(t==time[j]):
				l1[j,:2,i]=l1_error(f_lax[:N],time[j]),l1_error(f_upwind[:N],time[j])
				l1[j,2:,i]=l1_error(f_lax_wendroff[:N],time[j]),l1_error(f_lax_wendroff_v[:N],time[j])
				
for i in range(len(time)):
    plot_error(l1[i],i)
