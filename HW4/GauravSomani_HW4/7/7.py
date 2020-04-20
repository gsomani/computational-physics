import numpy as np
import math
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})
plt.rcParams["figure.figsize"] = (18,18)

d=0.02
alpha=math.pi/d

def T(u,w):
    if(u<=w):
    	return abs(math.sin(alpha*u))
    return 0

w=0.2
W=10*w
b=0.25
N=2**10
u=(np.arange(N)/N-0.5)*W
y=np.empty(N)

for i in range(N):
    y[i]=T(u[i],w)

c=np.fft.rfft(y)
xmax=50
kmax=round(xmax/b)
x=np.empty(1+2*kmax)
I=np.empty(1+2*kmax)

for i in range(1+kmax):
    I[kmax+i]=I[kmax-i]=(W*abs(c[i])/N)**2
    x[kmax-i]= b*i
    x[kmax+i]=-x[kmax-i]

plt.plot(x,I)
plt.grid(True, which="both")
plt.xticks(np.linspace(-50,50,21))
plt.xlabel('x (in mm)')
plt.ylabel('I(x)')
plt.title('Diffraction Intensity')
plt.savefig('diffraction.png') 
plt.show()
