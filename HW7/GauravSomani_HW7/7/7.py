import random
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt

def energy(L,J,S):
    E=0
    for i in range(N):
        for j in range(N):
            E += S[i][j]*(S[i-1][j]+S[i][j-1])
    E *= -J
    return E

def init(L):
    S = np.empty([L,L],dtype=int)
    for i in range(L):
        for j in range(L):
            S[i][j] = (-1)**random.randint(0,1)
    return S,np.sum(S)

def random_point(L):
    return random.randint(-1,L-2),random.randint(-1,L-2)

def metropolis(S,J,T,M,L):
    i,j = random_point(L)	
    diff = 2*J*S[i][j]*(S[i][j+1]+S[i][j-1]+S[i-1][j]+S[i+1][j])
    if(diff <= 0):
        S[i][j] *= -1
        return M+2*S[i][j]
    if(random.random()<=np.exp(-diff/T)):
        S[i][j] *= -1
        return M+2*S[i][j]
    return M

def monte_carlo(L,T,N):
    M = np.empty(N+1)
    S,M[0] = init(L)
    for i in range(N):
        M[i+1] = metropolis(S,1,T,M[i],L)
    return S,M

L=20
T=1
N=2**20

cmap = colors.ListedColormap(['red', 'blue'])
boundaries = [-1, 0, 1]
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

for T in range(1,4):
    S,M = monte_carlo(L,T,N)
    plt.plot(M)
    plt.title('Magnetisation at T = %i' % T)
    plt.ylabel('Magnetisation')
    plt.xlabel('time')
    plt.grid(True, which="both")
    plt.savefig('magnetisation_T'+str(T)+'.png')
    plt.show()
    plt.imshow(S,cmap=cmap,norm=norm)
    plt.colorbar()
    plt.title('T = %i' % T)
    plt.savefig('T = '+str(T)+'.png')
    plt.show()
