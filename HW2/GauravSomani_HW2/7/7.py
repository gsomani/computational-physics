import scipy as sp
import math
import matplotlib.pyplot as plt

gamma_unicode='\u0393'

def f(a,x):
    if(x==0):
        return 0
    return math.exp((a-1)*math.log(x)-x)

def gamma_diff(a,z):
    if(z==1):
        return 0
    x=z*(a-1)/(1-z)
    return (a-1)*f(a,x)/((z-1)*(z-1))

x=sp.linspace(0,5,101)

for a in range(2,5):
    g=[]
    for i in range(101):
        g.append(f(a,x[i]))
    plt.plot(x,g,label="a="+str(a))

plt.ylabel("$x^{a-1}e^{-x}$",fontsize=16)
plt.xlabel('x',fontsize=16)
plt.gca().legend()
plt.title("f(x)=$x^{a-1}e^{-x}$ for a=2,3,4")
plt.grid(True, which="both")
plt.savefig("fx.png")
plt.show()

def gamma(a):
    N=10000
    x=sp.linspace(0,1,N+1)
    total=0
    for i in range(1,N-1,2):
	    total+=4*gamma_diff(a,x[i])+2*gamma_diff(a,x[i+1])
    total+=gamma_diff(a,x[0])+4*gamma_diff(a,x[N-1])+gamma_diff(a,x[N])
    return total/(3*N)

print("At x=c, z=0.5  and c = a for peak of integrand to be at z=0.5")

print(gamma_unicode + "(1.5) = %.3f" % gamma(1.5) )

a=[3,6,10]

for i in a: 
    print(gamma_unicode + "(%i) = %.1f" % (i,gamma(i)) )
