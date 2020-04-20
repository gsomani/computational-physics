import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Circle

pi=constants.pi
epsilon=constants.epsilon_0
c=4*pi*epsilon
N=100
r=5
d=2*r
q=1
ratio=q/c

def potential_dipole(x,y):
    r1=((x-r)*(x-r)+y*y)**0.5
    r2=((x+r)*(x+r)+y*y)**0.5     
    return ratio*(r2-r1)/(r1*r2)

x=sp.linspace(-N//2,N//2,N+1)
p=sp.empty([N+1,N+1])
Exx=sp.empty([N+1,N+1])
Eyy=sp.empty([N+1,N+1])

for i in range(N+1):
    for j in range(N+1):
        p[i][j]=potential_dipole(x[j],x[i])

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

plt.text(-r,0,"-")
charge_colors = {True: '#aa0000', False: '#0000aa'}
plt.text(r,0,"+")
plt.title("Potential due to dipole")
plt.gca().set_aspect('equal')        
plt.pcolormesh(x, x, p,norm=colors.SymLogNorm(linthresh=ratio/(8*d),linscale=4))
plt.colorbar().set_label('Elecrostatic potential (V)')
plt.savefig("potential_dipole.png")
plt.show()

E = sp.log(sp.hypot(Exx, Eyy))
plt.streamplot(x, x, Exx, Eyy,color = E,linewidth=1, cmap=plt.cm.inferno,
              density=2, arrowstyle='->', arrowsize=1.5)

plt.gca().add_artist(Circle((-r,0), r/16, color=charge_colors[0]))
plt.gca().add_artist(Circle((r,0), r/16, color=charge_colors[1]))

plt.xlim(-N//2,N//2)
plt.ylim(-N//2,N//2)
plt.gca().set_aspect('equal')
plt.text(-r,0,"-")
plt.text(r,0,"+")
plt.colorbar().set_label('Elecric Field (V/cm)')
plt.title("Electric field due to dipole")
plt.savefig("Electric_field_dipole.png")
plt.show()
