import scipy as sp
import math
import matplotlib.pyplot as plt
from gaussxw import gaussxw,gaussxwab

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,15)

G=6.674e-11
sigma=100
M=500
N=100
L=10
x,w=gaussxwab(N,0,L//2)

def f(x,y,z):
    return (x*x+y*y+z*z)**(-1.5)

def Fz(z):
    total=0
    for i in range(N):
        for j in range(N):
            total+=w[i]*w[j]*f(x[i],x[j],z)
    return 4*G*sigma*z*total

F=sp.empty(M+1)
z=sp.linspace(0,L,M+1)
F[0]=2*sp.pi*sigma*G
for i in range(1,M+1):
    F[i]=Fz(z[i])

plt.plot(z[1:],F[1:])
plt.plot(z[0],F[0],'.',markerfacecolor='none',markersize=8)
plt.xticks(sp.linspace(0,10,21))
plt.ylabel("$F_z$(z) [in N]",fontsize=16)
plt.xlabel('z [in m]',fontsize=16)
plt.title("Force in z-direction on 1kg point mass")
plt.grid(True, which="both")
plt.savefig("Fz.png")
plt.show()
