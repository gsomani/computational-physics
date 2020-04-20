import random
import numpy as np
import matplotlib.pyplot as plt

def gaussian_random(sigma):
	r=(-2*sigma*sigma*np.log(1-random.random()))**0.5
	return r*np.cos(np.pi*random.random())

def f(x):
	return x*x - np.cos(4*np.pi*x)

def metropolis(x,b,sigma,E):
	d = gaussian_random(sigma)	
	diff = E(x+d) - E(x)
	if(diff < 0):
		return x+d
	if(np.log(random.random())<=-b*diff):
		return x+d
	return x

def simulated_annealing(xi,bi,bf,tau,sigma,E):
    b,x=bi,xi
    X = [xi]
    while(b<=bf):	
        x = metropolis(x,b,sigma,E)
        X.append(x)
        b += b/tau
    fig, ax = plt.subplots()
    plt.plot(X,".",markersize=1)
    plt.xlabel('time')
    plt.ylabel('x')
    plt.grid(True, which="both")
    return x,[fig,ax]

infinity = 2**30
cooling_time = 2**14
bi = 2**-2

xi = 2
x,[fig,ax] = simulated_annealing(xi,bi,infinity,cooling_time,1,f)
plt.title('Simulated annealing { $f(x) = x^2 − cos(4 \pi x)$ }')
plt.savefig('a.png')
plt.show()
print("a) x (by simulated annealing) = %.1E" %x)

root_2 = 2**0.5
root_3 = 3**0.5

def g(x):
    if(x<=0 or x>=50):
        return 2**1023
    return np.cos(x) + np.cos(root_2*x) + np.cos(root_3*x)

xr = [0,50]
xi = 0.5*(xr[0]+xr[1])
x,[fig,ax] = simulated_annealing(xi,bi,infinity,cooling_time,10,g)
plt.title('Simulated annealing { $f(x) = cos(x) + cos(\sqrt{2} x) + cos(\sqrt{3} x)$ }')
plt.savefig('b.png')
plt.show()
print("b) x (by simulated annealing) = %.2f" %x)
