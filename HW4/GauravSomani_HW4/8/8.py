import numpy as np
import matplotlib.pyplot as plt
import math
from cmath import phase
from numpy.fft import rfft2,irfft2,fft2

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,15)

file = open('blur.txt', 'r')
n=1024
b=np.empty([n,n])

i=0
for line in file:
    b[i]=line.split()
    i+=1

plt.imshow(b,'gray')
plt.title('Blurred image')
plt.savefig('blur.png')
plt.show()

sigma=25
root_2pi=(2*math.pi)**0.5

c=np.empty([n,n])

def fs(x,y,sigma):
    return math.exp(-(x*x+y*y)/(2*sigma*sigma))/(sigma*root_2pi)

f=np.empty([n,n])

for i in range(1+n//2):
    for j in range(1+n//2):
        f[i][j]=f[-i][j]=f[i][-j]=f[-i][-j]=fs(i,j,sigma)

plt.imshow(f,'gray')
plt.title('Point spread function')
plt.savefig('gaussian.png')
plt.show()

B=rfft2(b)
F=fft2(f)

A=np.empty([n,1+n//2],complex)
eps=2**-10
p=[]

for i in range(n):
    for j in range(1+n//2):
        if(abs(F[i][j])>eps):
            A[i][j]=B[i][j]/F[i][j]
        else:
            A[i][j]=B[i][j]

a=irfft2(A)
plt.imshow(a,'gray')
plt.title('Deblurred image')
plt.savefig('deblur.png')
plt.show()
