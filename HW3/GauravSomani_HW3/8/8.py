import matplotlib.pyplot as plt
import scipy as sp

def P(a,x):
    l=len(a)
    y=a[0]
    for i in range(1,l):
        y+=a[i]*x**i
    return y

def deriv_P(a):
    l=len(a)
    y=[a[1]]
    for i in range(2,l):
        y.append(i*a[i])
    return y

def newton(guess,a,b,eps):
    f=P(a,guess)
    deriv_f=P(b,guess)
    x=guess-f/deriv_f
    g1=abs(P(a,x+0.5*eps))
    g2=abs(P(a,x-0.5*eps))
    g=abs(P(a,x))
    if(min(g1,g2)>g):
        return x
    return newton(x,a,b,eps)    

a=[1,-42,420,-1680,3150,-2772,924]
b=deriv_P(a)
x=sp.linspace(0,1,201)
y=sp.empty(201)

for i in range(201):
    y[i]=P(a,x[i])

plt.plot(x,y)
plt.grid(True,"both")
plt.title("P(x) = 924$x^6$ - 2772$x^5$ + 3150$x^4$ - 1680$x^3$ + 420$x^2$ - 42$x$ + 1")
plt.savefig("polynomial.png")
plt.show()

root=[0.05,0.15,0.35,0.65,0.85,0.95]
eps=1e-10

for i in range(6):
    root[i]=round(newton(root[i],a,b,eps),11)

print("Roots = "+str(root))
