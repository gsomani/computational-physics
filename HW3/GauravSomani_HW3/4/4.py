import scipy as sp
import math
import matplotlib.pyplot as plt

def f(c,x):
    return 1-math.exp(-c*x)

def fixed_point(guess,c,eps,p):
    global it
    x=f(c,guess)
    it+=1    
    g1=abs(f(c,x+0.5*eps)-(x+0.5*eps))
    g2=abs(f(c,x-0.5*eps)-(x-0.5*eps))
    g=abs(f(c,x)-x)
    if(p=='print'):
        print("iteration %i, x = %.6f" % (it,x))    
    if(min(g1,g2)>g):
        return x
    return fixed_point(x,c,eps,p)

def solve_x(c,eps):
    if(c==0):
        return 0
    guess=max(0,2*(c-1)/(c*c))
    return fixed_point(guess,c,eps,0)

it=0
c=2
x0=1
eps=1e-6
print("Fixed point iteration without acceleration")
print("\nInitial guess, x = %.6f" % x0)
xc=fixed_point(x0,c,eps,'print')
print("\nx = %.6f , iterations = %i" %(xc,it))

N=300
c=sp.linspace(0,3,N+1)
x=sp.empty(N+1)
    
for i in range(N+1):
    x[i]=solve_x(c[i],eps)

plt.plot(c,x)
plt.xlabel('c',fontsize=16)
plt.ylabel('x',fontsize=16) 
plt.title("Percolation Transition")       
plt.yticks(sp.linspace(0,1,11))
plt.xticks(sp.linspace(0,3,16))
plt.grid(True,"both")
plt.savefig("percolation.png")
plt.show()

def g(r,c,x):
    return x+r*(f(c,x)-x) 

def fixed_point_accelerated(r,guess,c,eps):
    global it
    x=g(r,c,guess)
    it+=1
    g1=abs(f(c,x+0.5*eps)-(x+0.5*eps))
    g2=abs(f(c,x-0.5*eps)-(x-0.5*eps))
    gx=abs(f(c,x)-x)
    print("iteration %i, x = %.6f" % (it,guess))    
    if(min(g1,g2)>gx):
        return x
    return fixed_point_accelerated(r,x,c,eps)

it=0
c=2
x0=1
print("\n\nFixed point iteration with acceleration")
print("\nInitial guess, x = %.6f" % x0)
xc=fixed_point_accelerated(1.5,x0,c,eps)
print("x = %.6f , iterations = %i" %(xc,it))
