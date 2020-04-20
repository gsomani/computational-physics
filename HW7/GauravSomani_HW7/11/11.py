import random
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

def random_dir():
    r = random.randint(0,3)
    a,b = r//2 , 2*(r%2)-1
    if(a):
        return 0,b
    return b,0

def brownian_update(x,y,L,P):
    diff = random_dir()
    x_,y_ = x+diff[0],y+diff[1]
    stuck = False
    if(x_==0 or x_==L-1 or y_==0 or y_==L-1):
        stuck = True
    elif(P[x_+1][y_]>=1 or P[x_-1][y_]>=1 or P[x_][y_+1]>=1 or P[x_][y_-1]>=1):
        stuck = True
    return [x_,y_],stuck

def DLA_1(L):
    P = np.zeros([L,L],dtype=int)
    x = [0,0]
    center = [L//2,L//2]
    i=j=0
    while(x!=center):
        stuck=False
        x = center
        S = np.copy(P)
        while(stuck==False):
            x,stuck = brownian_update(x[0],x[1],L,P)
        i+=1
        P[x[0]][x[1]]=i
        if(i%(2**j)==0):
            plt.imshow(P)
            plt.title('Diffusion-limited aggregation at edges {Number of particles = %i}' % i)
            plt.savefig('edge'+str("%02d" % j)+'.png')
            j+=1
    plt.imshow(P)
    plt.title('Diffusion-limited aggregation at edges {Number of particles = %i}' % i)
    plt.savefig('edge'+str("%02d" % j)+'.png')
    return P

def random_point(c,d):
    x,y = c[0]+d*np.cos(2*np.pi*random.random()),c[1]+d*np.sin(2*np.pi*random.random())
    return [int(np.ceil(x)),int(np.ceil(y))]

def DLA_2(L):
    P = np.zeros([L,L],dtype=int)
    x = [0,0]
    center = [L//2,L//2]
    r = 0
    P[center[0]][center[1]]=1
    i=j=0
    while(2*r<L//2):
        stuck=False
        x = random_point(center,r+1)
        while(stuck==False):
            x,stuck = brownian_update(x[0],x[1],L,P)
            dist = np.ceil(la.norm([x[0]-center[0],x[1]-center[1]]))
            if(dist>=2*(r+1)):
                x = random_point(center,r+1)
                stuck = False
        i+=1
        P[x[0]][x[1]]=i
        r = max(int(dist),r)
        if(i%(2**j)==0):
            plt.imshow(P)
            plt.title('Diffusion-limited aggregation at center {Number of particles = %i}' % i)
            plt.savefig('center'+str("%02d" % j) +'.png')
            j+=1
    plt.imshow(P)
    plt.title('Diffusion-limited aggregation at center {Number of particles = %i}' % i)
    plt.savefig('center'+str("%02d" % j)+'.png')
    return P

L=101

P = DLA_1(L)
plt.imshow(P)
plt.title('Diffusion-limited aggregation at edges')
plt.savefig('edge_dla.png')
plt.show()

P = DLA_2(L)
plt.imshow(P)
plt.title('Diffusion-limited aggregation at center')
plt.savefig('center_dla.png')
plt.show()
