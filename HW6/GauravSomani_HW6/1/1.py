import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

N=100
L=10
r=[2,8]
V=1

def set_region(reg,f,N,L,r,V):
    g=[r[0]*N//L,r[1]*N//L]
    for i in range(g[0],g[1]+1):
        reg[i][g[0]]=1
        f[i][g[0]]=V
        reg[i][g[1]]=-1
        f[i][g[1]]=-V

def relaxation(reg,f,N,w):
    err=0
    for i in range(1,N):
        for j in range(1,N):
            if(reg[i][j]==0):
                diff = w*((f[i-1][j]+f[i+1][j]+f[i][j-1]+f[i][j+1])/4 - f[i][j])
                f[i][j]+=diff
            err=max(abs(diff),err)
    return f,err

eps=2**-20

def solve(eps,N,L,V):
    reg = np.zeros([N+1,N+1])
    f = np.zeros([N+1,N+1])
    set_region(reg,f,N,L,r,V)
    err=eps
    w=2/(1+np.sin(np.pi/N))
    while(err>eps/4):
        f,err=relaxation(reg,f,N,w)
    return f,np.linspace(0,L,N+1)
    
phi,x = solve(eps,N,L,V)
plt.pcolormesh(x,x,phi)
plt.xlabel('width (in cm)')
plt.ylabel('height (in cm)')
plt.gca().set_aspect('equal')
plt.colorbar().set_label('Electrostatic potential (V)')
plt.title('Electrostatic potential inside a capacitor')
plt.savefig('potential.png')
plt.show()
