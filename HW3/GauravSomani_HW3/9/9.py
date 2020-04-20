import matplotlib.pyplot as plt
import scipy as sp

m=7.348 # in 10^22 kg
M=597.4 # in 10^22 kg
b=m/M
i=0

def P(a,x):
    l=len(a)
    y=a[0]
    for i in range(1,l):
        y+=a[i]*x**i
    return y

def secant(guess,a,eps):
    x0=guess[0]
    x1=guess[1]    
    f0=P(a,x0)
    f1=P(a,x1)
    global i
    i+=1
    x=(x1*f0-x0*f1)/(f0-f1)
    g1=abs(P(a,x+0.5*eps))
    g2=abs(P(a,x-0.5*eps))
    g=abs(P(a,x))
    if(min(g1,g2)>g):
        return x
    return secant([x1,x],a,eps)    

a=[-1,2,b-1,1,-2,1]
guess=[0.8,0.9]
eps=1e-5

R=3.844e8
r=secant(guess,a,eps)*R
print("r = %.3E" % r)

