import scipy as sp
from scipy import constants
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Circle

pi=constants.pi
epsilon=constants.epsilon_0
c=4*pi*epsilon
N=100
q0=1e-2
L=10
k=2*pi/L
x=sp.linspace(-N//2,N//2,N+1)
p=sp.empty([N+1,N+1])
Exx=sp.empty([N+1,N+1])
Eyy=sp.empty([N+1,N+1])
M=20
xx=sp.linspace(-L//2,L//2,M+1)

def f(x,y,x0,y0):
    r=((x-x0)*(x-x0)+(y-y0)*(y-y0))**0.5
    if(2*r<L/M):
        return 0
    return q0*math.sin(k*x)*math.sin(k*y)/r

def w(i):
    if(i==0 or i==M):
        return 0.5
    return 1

def Vz(x,y):
    total=0
    for i in range(M+1):
        for j in range(M+1):
            total+=w(i)*w(j)*f(xx[j],xx[i],x,y)
    return total/(M*c)

for i in range(N+1):
    for j in range(N+1):
        p[i][j]=Vz(x[j],x[i])

print(p[N//2][N//2])
def E(i,j):
    if(i==0):
        Ey=p[i][j]-p[i+1][j]
    else: 
        if(i==N): 
            Ey=p[i-1][j]-p[i][j]
        else:
            Ey=(p[i-1][j]-p[i+1][j])/2
    if(j==0):
        Ex=p[i][j+1]-p[i][j]
    else: 
        if(j==N): 
            Ex=p[i][j-1]-p[i][j]
        else:
            Ex=(p[i][j-1]-p[i][j+1])/2
    return Ex,Ey
    
for i in range(N+1):
    for j in range(N+1):
        Exx[i][j],Eyy[i][j]=E(i,j)

# Plot the streamlines with an appropriate colormap and arrow style
color = sp.hypot(Exx, Eyy)
plt.streamplot(x, x, Exx, Eyy, color=color,linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='->', arrowsize=1.5)

plt.xlim(-N//2,N//2)
plt.ylim(-N//2,N//2)
plt.gca().set_aspect('equal')
plt.colorbar().set_label('Elecric Field (V/cm)')
plt.title("Electric field due to $\sigma(x,y)$=$q_0 sin(2\pi x/L) sin(2\pi y/L)$")
plt.savefig("Electric_field_sine.png")
plt.show()
