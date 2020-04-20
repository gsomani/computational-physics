import scipy as sp
from scipy import constants
import math
import cmath
import matplotlib.pyplot as plt

pi=constants.pi
d=20
alpha=pi/d
lamda=0.5
f=1e6
j=complex(0,1)

def diff(u,q,x):
    return cmath.exp(j*2*pi*x*u/(lamda*f))*(q**0.5)

def int1d(q,m):
    n=500
    x=sp.linspace(-5e4,5e4,n+1)
    intensity=sp.empty(n+1)
    du=10*d/m
    for k in range(n+1):
        total=0
        for i in range(1,m-1,2):    
            total+=4*diff(u[i],q[i],x[k])+2*diff(u[i+1],q[i+1],x[k])  
        total+=4*diff(u[m-1],q[m-1],x[k])+diff(u[0],q[0],x[k])+diff(u[m],q[m],x[k])
        total*=du/3
        intensity[k]=abs(total)*abs(total)
    pattern=sp.empty([100,n+1])
    for i in range(100):
        pattern[i]=intensity
    return pattern

m=200
u=sp.linspace(-5*d,5*d,m+1)    
q=sp.empty(m+1)

def q1(u):
	return math.sin(alpha*u)**2

for i in range(m+1):
    q[i]=q1(u[i])    
plt.imshow(int1d(q,m),'gray')
plt.title('q(u)=$sin^2$($\\alpha$u)')
plt.colorbar()
plt.savefig('diffraction_1.png')
plt.show()

def q2(u):
	return math.sin(alpha*u)**2 * math.sin(0.5*alpha*u)**2

for i in range(m+1):
    q[i]=q2(u[i])
plt.imshow(int1d(q,m),'gray')
plt.title('q(u)=$sin^2$($\\alpha$u)$sin^2$($\\beta$u)')
plt.colorbar()
plt.savefig('diffraction_2.png')
plt.show()

def q3(u):
    if(u<=-35 or u>=25):
        return 1
    return 0	

m=180
u=sp.linspace(-45,45,m+1)    
q=sp.empty(m+1)

for i in range(m+1):
    q[i]=q3(u[i])
plt.imshow(int1d(q,m),'gray')
plt.colorbar()
plt.title('Square Slits')
plt.savefig('diffraction_3.png')
plt.show()
