import scipy as sp
from scipy import constants
import math
import matplotlib.pyplot as plt
from gaussxw import gaussxw,gaussxwab

pi=constants.pi
sqrt_pi=pi**0.5

def H(n,x):
	h=[1,2*x]	
	for i in range(2,n+1):
		h.append(2*x*h[i-1]-2*(i-1)*h[i-2])
	return h[n]

def psi(n,x):
	return math.exp(-x*x/2)*H(n,x)/((math.factorial(n)*sqrt_pi*(1<<n))**0.5)

x=sp.linspace(-4,4,161)
ps=[]
for j in range(4):
	ps=[]
	for i in range(161):
		ps.append(psi(j,x[i]))
	plt.plot(x,ps,label="$\psi_"+str(j)+"$")

plt.ylabel("$\psi$(x)",fontsize=16)
plt.xlabel('x',fontsize=16)        
plt.gca().legend()
plt.title("$\psi_n$(x) for n=0,1,2,3")
plt.grid(True, which="both")
plt.savefig("wavefunctions.png")
plt.show()

N=800
x=sp.linspace(-10,10,N+1)
ps=sp.empty(N+1)
for i in range(N//2,N+1):
	ps[i]=ps[N-i]=psi(30,x[i])

plt.ylabel("$\psi$(x)",fontsize=16)
plt.xlabel('x',fontsize=16)        
plt.plot(x,ps)
plt.title("$\psi_{30}$(x)")
plt.grid(True, which="both")
plt.savefig("30th_state.png")
plt.show()

def expectation_function(n,x):
    return ((x*H(n,x))**2)*math.exp(-x*x)/(math.factorial(n)*sqrt_pi*(1<<n))

def g(n,z):
    c=(n+1)**0.5
    if(z==1):
        return 0  
    x=c*z/(1-z)
    return c*expectation_function(n,x)/((1-z)*(1-z))

def uncertainity(n):
    x,w=gaussxwab(100,0,1)
    total=0
    for i in range(100):
        total+=w[i]*g(n,x[i])
    return (2*total)**0.5

print("Uncertainity for n=5 is %.6f" % uncertainity(5))
