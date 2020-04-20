import scipy as sp
from scipy import constants
import matplotlib.pyplot as plt
import math

pi=constants.pi
points=5000
x=sp.zeros(points)
y=sp.zeros(points)
plt.grid(True, which="both")

def parametric_xy(theta):
    x=2*math.cos(theta)+math.cos(2*theta)
    y=2*math.sin(theta)-math.sin(2*theta)
    return x,y

theta=sp.linspace(0,2*pi,points)

for i in range(points):
    x[i],y[i]=parametric_xy(theta[i])
plt.plot(x,y)
plt.title('(a)')
plt.savefig('a.png')
plt.show()

def polar_xy(r,theta):
    x=r*math.cos(theta)
    y=r*math.sin(theta)
    return x,y

theta=sp.linspace(0,10*pi,points)
for i in range(points):
    r=theta[i]*theta[i]
    x[i],y[i]=polar_xy(r,theta[i])

plt.grid(True, which="both")
plt.plot(x,y)
plt.title('(b)')
plt.savefig('b.png')
plt.show()

theta=sp.linspace(0,24*pi,points)
for i in range(points):
    r=math.exp(math.cos(theta[i]))-2*math.cos(4*theta[i])+math.sin(theta[i]/12)**5
    x[i],y[i]=polar_xy(r,theta[i])

plt.grid(True, which="both")
plt.title('(c)')
plt.plot(x,y)
plt.savefig('c.png')
plt.show()
