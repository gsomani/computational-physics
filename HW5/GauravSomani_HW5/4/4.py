import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 24})
plt.rcParams["figure.figsize"] = (25,25)

def rk4_step(x,t,h,f):
	k0=f(x,t)
	k1=f(x+0.5*h*k0,t+0.5*h)
	k2=f(x+0.5*h*k1,t+0.5*h)
	k3=f(x+h*k2,t+h)
	xf = x+h*(k0+2*k1+2*k2+k3)/6
	return xf,t+h

def rk4(x0,time,h,f):
    d=len(x0)
    N=int(round(1+time/h)) 
    x,t=np.empty([N,d]),np.empty(N)
    x[0],t[0]=x0,0
    for i in range(1,N):
    	x[i],t[i]=rk4_step(x[i-1],t[i-1],h,f)
    return x,t

w=1

def f_shm(R,t):
    fx=R[1]
    fv=-w*w*R[0]
    return np.array([fx,fv])  

def f_anharmonic(R,t):
    fx=R[1]
    fv=-w*w*(R[0]**3)
    return np.array([fx,fv])  

def f_pol(R,t):
    fx=R[1]
    fv = u*(1-R[0]*R[0])*R[1]-w*w*R[0]
    return np.array([fx,fv])  

def phase_space_plot(x,v,l,lab,title):
    plt.plot(x[:l],v[:l],label=lab)
    plt.grid(True, which="both")
    plt.title('Phase Space ($\omega = 1$)'+title)
    plt.xlabel('x')
    plt.ylabel('v')

def position_plot(t,x,lab,title):
    plt.plot(t,x,label=lab)
    plt.grid(True, which="both")
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.title('Trajectory ($\omega = 1$)'+title)

f=[f_shm,f_anharmonic]
amp = [2,1]
filename=['shm_trajectory.png','shm_phase.png','anharmonic_trajectory.png','anharmonic_phase.png']
title=['{Simple Harmonic Oscillator}','{Anharmonic Oscillator}']

for j in range(2):
	for i in range(2):
		R,t=rk4([amp[i],0],50,2**-10,f[j])
		x,v = R[:,0],R[:,1]
		position_plot(t,x,'Amplitude = '+str(amp[i]),title[j])
	plt.gca().legend()
	plt.savefig(filename[2*j])
	plt.show()    
	phase_space_plot(x,v,2**13,'',title[j])
	plt.savefig(filename[2*j+1])
	plt.show()
    
mu=[1,2,4]
for i in range(3):
    u=mu[i]
    R,t=rk4([1,0],20,2**-10,f_pol)
    x,v=R[:,0],R[:,1]
    phase_space_plot(x,v,len(x),'$\mu = '+str(u)+'$','{van der Pol oscillator}')

plt.gca().legend()
plt.savefig('van_der_pol.png')
plt.show()
