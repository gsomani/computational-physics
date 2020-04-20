import scipy as sp
import matplotlib.pyplot as plt

approx= '\u2248'
plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,15)

step=0.01
points=3/step
r=sp.linspace(1,4,points)
x=sp.linspace(0.5,0.5,points)
for i in range(1000):
    x=r*x*(1-x)

plt.ylabel('x',fontsize=16)
plt.xlabel('r',fontsize=16)    
plt.xticks(sp.linspace(1, 4, 16))
plt.yticks(sp.linspace(0, 1, 21))
plt.grid(True, which="both")


for i in range(1000):
    x=r*x*(1-x)
    plt.plot(r,x,'k.',markersize=1)

print("(a)")
print("For a fixed point x, x  = rx(1-x). So, fixed point lies on curve x=(r-1)/r. For a fixed point, x takes a single value (a single point on the plot on vertical line at r) for that r. For limit cycle, x takes a finite number of values (usually 2,4 or 8 values) denoted by finite number of points on vertical line at r. For chaos at a particular r, there will be large number of points (number of points should always increase forming pseudo-continuous vertical line) on the vertical line at r.\n")  
print("(b)")
print("r "+ approx + " 3.57 is the edge of the chaos.")
plt.annotate('edge of chaos', xy=(3.57, 0.35), xytext=(3.57,-0.05),arrowprops=dict(facecolor='black', shrink=0.05))
plt.annotate('fixed point', xy=(2, 0.5), xytext=(2,0.3),arrowprops=dict(facecolor='black', shrink=0.05))
plt.annotate('limit cycle', xy=(3.2, 0.8), xytext=(3.2,0.9),arrowprops=dict(facecolor='black', shrink=0.05))
plt.savefig('figtree.png')    
plt.show()
