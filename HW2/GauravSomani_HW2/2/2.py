import scipy as sp
from scipy import constants
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors

pi=constants.pi

def J(m,x):
	N=1000
	theta=sp.linspace(0,pi,N+1)
	total=0
	for i in range(1,N-1,2):
		total+=4*math.cos(m*theta[i]-x*math.sin(theta[i]))+2*math.cos(m*theta[i+1]-x*math.sin(theta[i+1]))
	total+=math.cos(m*theta[0]-x*math.sin(theta[0]))+4*math.cos(m*theta[N-1]-x*math.sin(theta[N-1]))+math.cos(m*theta[N]-x*math.sin(theta[N]))
	return total/(3*N)

def I(r,k):
	if(r==0):
		return 0.25
	return (J(1,k*r)/(k*r))**2

x=sp.linspace(0,20,101)
J0=sp.linspace(0,20,101)
J1=sp.linspace(0,20,101)
J2=sp.linspace(0,20,101)

for i in range(101):
    J0[i]=J(0,x[i])
    J1[i]=J(1,x[i])
    J2[i]=J(2,x[i])

plt.ylabel('J(x)',fontsize=16)
plt.xlabel('x',fontsize=16)        
plt.plot(x,J0,label="$J_0$")
plt.plot(x,J1,label="$J_1$")
plt.plot(x,J2,label="$J_2$")
plt.grid(True, which="both")
plt.gca().legend()
plt.title('Bessel functions')
plt.savefig('bessel.png')
plt.show()

lamda=0.5
k=2*pi/lamda
N=200
x=sp.linspace(-1,1,N+1)
y=sp.linspace(-1,1,N+1)
Intensity=sp.empty([N+1, N+1])

for i in range(N//2,N+1):
	for j in range(N//2,N+1):
		r=(x[j]*x[j]+y[i]*y[i])**0.5
		Intensity[N-i][N-j]=Intensity[N-i][j]=Intensity[i][N-j]=Intensity[i][j]=I(r,k)

title = ['Diffraction Pattern (Intensity)','Diffraction Pattern (Intensity on Logarithm scale)']
filename = ['diffraction_intensity.png','log_intensity.png']
norm_type = [ None, colors.LogNorm()]

for i in range(2):
    plt.gca().set_aspect('equal')
    plt.pcolormesh(x,y,Intensity,norm = norm_type[i] )
    plt.title(title[i])
    plt.colorbar().set_label('Intensity')
    plt.savefig(filename[i])
    plt.show()
