import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

theta_unicode='\u03B8'
xi_unicode='\u03BE'
h=2**-10
eps=h/2
n=3
r0=2

def rk4_step(g,x,h,n,f):
	k0=f(g,x,n)
	k1=f(g+0.5*h*k0,x+0.5*h,n)
	k2=f(g+0.5*h*k1,x+0.5*h,n)
	k3=f(g+h*k2,x+h,n)
	diff = h*(k0+2*k1+2*k2+k3)/6
	return diff

def rk4(g0,x0,h,f,n,eps):
    g=g0
    x=x0
    m,theta=[g0[0]],[g0[1]]
    X=[x0]
    h0=h
    while(g[1]>eps):
        diff = rk4_step(g,x,h,n,f)
        if(diff[1] < -g[1]):
        	h=h/2
        	continue
        g+=diff 	
        m.append(g[0])
        theta.append(g[1])	
        x+=h
        X.append(x)
    radius=x
    return m,theta,X,radius

def f(r,x,n):
    fm=x*x*r[1]**n
    f_theta=-r[0]/(x*x)
    return np.array([fm,f_theta])

def secant_root(x,f,eps):
	while(abs(x[1]-x[0])>=eps/2):
		y=[f(x[0]),f(x[1])]
		m=(y[1]-y[0])/(x[1]-x[0])
		x = x[1],x[1]-y[1]/m
	return x[1]

def diff_rad(theta):
	g0=[(eps*theta)**3/3,theta]
	m,th,x,radius=rk4(g0,eps,eps,f,n,eps)
	return radius-r0

M=19
rad_surface=np.empty(M)
theta_core=np.linspace(1,10,M)
j=0
for i in range(0,M):
    theta_0 = theta_core[i]
    g0=[(eps*theta_0)**3/3,theta_0]
    m,theta,x,radius=rk4(g0,eps,h,f,n,eps)
    rad_surface[i] = radius 
    if(radius > r0):
        j=i
    if(i%2==0):
        plt.plot(x,theta,label='$\\theta_0=$'+str(int(theta_0)))
plt.grid(True,"both")	
plt.xlabel('$\\xi$')
plt.ylabel('$\\theta$')
plt.gca().legend()
plt.yticks(np.linspace(0,10,11))
plt.xticks(np.linspace(0,10,11))
plt.title('Solution to Lane-Emden equations')
plt.savefig('shooting_hand.png')
plt.show()

plt.plot(theta_core,rad_surface,label='$\\theta_0=$'+str(theta_0))
plt.grid(True,"both")	
plt.ylabel('$\\xi^*$')
plt.xlabel('$\\theta_0$')
plt.xticks(np.linspace(1,10,19))
plt.title('$\\xi^*$ vs $\\theta_0$')
plt.savefig('density.png')
plt.show()

print('From plot,'+theta_unicode+' is between %.1f and %.1f'%(theta_core[j],theta_core[j+1]))
theta_0=secant_root([theta_core[j],theta_core[j+1]],diff_rad,eps)
print('Required '+theta_unicode+'_0 = %.3f for ' % theta_0 + xi_unicode+'* = %i '%r0)

g0=[(eps*theta_0)**3/3,theta_0]
m,theta,x,radius=rk4(g0,eps,h,f,n,eps)
plt.grid(True,"both")	
plt.xlabel('$\\xi$')
plt.ylabel('$\\theta$')
plt.plot(x,theta)
plt.title('Solution corresponding to $\\xi^* = 2$')
plt.savefig('solution.png')
plt.show()
