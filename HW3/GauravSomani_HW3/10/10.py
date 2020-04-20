import math

I0=3e-6
R1=1
R2=4
R3=3
R4=2
Vt=0.05
V=5

def f(x,y):
    return (1/R1+1/R2)*x+I0*(math.exp((x-y)/Vt)-1)-V/R1

def g(x,y):
    return (1/R3+1/R4)*y-I0*(math.exp((x-y)/Vt)-1)-V/R3

def fy(x,y):
    return -I0*math.exp((x-y)/Vt)/Vt

def gx(x,y):
    return fy(x,y)

def fx(x,y):
    return (1/R1+1/R2)-fy(x,y)

def gy(x,y):
    return (1/R3+1/R4)-fy(x,y)

def newton_iteration(x,y):
    det=fx(x,y)*gy(x,y)-gx(x,y)*fy(x,y)
    diff_x=(fy(x,y)*g(x,y)-gy(x,y)*f(x,y))/det
    diff_y=(gx(x,y)*f(x,y)-fx(x,y)*g(x,y))/det
    return x+diff_x,y+diff_y

def newton(guess,eps):
    x,y=newton_iteration(guess[0],guess[1])
    diff=[x-guess[0],y-guess[1]]
    diff_size=diff[0]*diff[0]+diff[1]*diff[1]
    if(diff_size<0.25*eps*eps):
        return x,y
    return newton([x,y],eps)

Vi=[V*R2/(R1+R2),V*R3/(R3+R4)] 
eps=1e-6

Vi=newton(Vi,eps)

print("V1 = %.3f V ,V2 = %.3f V" %(Vi[0],Vi[1]))
