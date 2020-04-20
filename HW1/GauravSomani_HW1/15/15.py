import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math

def outside_mandelbrot(c):
    z=complex(0,0)
    for i in range(100):
        z=z*z+c
        if(abs(z)>2):
            return i+1,True
    return 100,False

N=500
x=sp.linspace(-2,2,N+1)
p=sp.empty([N+1, N+1]) 

it=sp.empty([N+1, N+1])

for i in range(N+1):
    for j in range(N+1):
        c=complex(x[j],x[i])
        it[i][j],p[i][j]=outside_mandelbrot(c)

plt.gca().set_aspect('equal')
plt.pcolormesh(x, x, p,cmap="gray")
plt.title('Mandelbrot set')
plt.savefig('mandelbrot.png')
plt.show()
plt.gca().set_aspect('equal')
plt.pcolormesh(x, x, it,cmap="jet")
plt.title('Mandelbrot set (Number of iterations)')
plt.colorbar().set_label('Number of iterations')
plt.savefig('mandelbrot_iterations.png')
plt.show()
plt.gca().set_aspect('equal')        
plt.pcolormesh(x, x, it,cmap="jet",norm=colors.LogNorm())
plt.title('Mandelbrot set (Number of iterations on logarithm scale)')
plt.colorbar().set_label('Number of iterations')
plt.savefig('mandelbrot_log_iterations.png')
plt.show()
