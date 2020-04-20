import random
import numpy as np
import matplotlib.pyplot as plt

def random_bin(L):
    p = random.randint(0,L*L-1)
    x,y = p//L , p%L
    m = [x,y]
    n = [-1,-1]
    while(n[0]<0 or n[0]>=L or n[1]<0 or n[1]>=L):
        r = random.randint(0,3)
        a,b = r//2 , 2*(r%2)-1
        if(a):
            dx,dy = 0,b
        else:
            dx,dy = b,0
        n = [x+dx,y+dy]
    return m,n

def remove_dimer(index,lat,dim):
    N = len(dim)
    p,q = dim[index],dim[N-1]
    lat[q[0][0]][q[0][1]] = lat[q[1][0]][q[1][1]] = index 
    lat[p[0][0]][p[0][1]] = lat[p[1][0]][p[1][1]] = -1
    dim[index],dim[N-1] = q,p
    dim.pop()
    return

def add_dimer(p,q,lat,dim):
    lat[p[0]][p[1]] = lat[q[0]][q[1]] = len(dim)
    dim.append([p,q])
    return

def metropolis(b,lat,dim):
    p,q = random_bin(L)
    r,s = lat[p[0]][p[1]],lat[q[0]][q[1]]
    diff = r - s
    dE=0
    if(diff==0):
        if(r == -1 and s == -1):
            add_dimer(p,q,lat,dim)
            dE -= 1
        else:
            if(np.log(random.random()) <= -b):
                remove_dimer(r,lat,dim)
                dE += 1 
    return dE

def simulated_annealing(bi,bf,L,tau):
    lat = -1*np.ones([L,L],dtype=int)
    dim = []
    b = bi
    E=i=0
    j=0
    while(b<=bf):	
        dE = metropolis(b,lat,dim)
        E += dE
        b += b/tau
        i+=1
        if(i==(2**j)):
            plt.imshow(lat>=0,cmap='Greys')
            plt.title('Dimer Covering {cooling time constant = '+str(tau)+' } (t = ' + str(i) + ')')
            plt.savefig(str("%i_%02d" % (tau,j))+'.png')
            j+=1
    plt.imshow(lat>=0,cmap='Greys')
    plt.title('Dimer Covering {cooling time constant = '+str(tau)+' } (t = ' + str(i) + ')')
    plt.savefig(str("%i_%02d" % (tau,j))+'.png')
    return E

infinity = 750
cooling_time = [2**7,2**10,2**13,2**16]
bi = 2**-2
L = 50
for i in range(4):
    x = -simulated_annealing(bi,infinity,L,cooling_time[i])
    print("Number of dimers (cooling time constant = %i) = %i" %(cooling_time[i],x)) 
