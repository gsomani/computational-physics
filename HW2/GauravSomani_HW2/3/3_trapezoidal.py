import scipy as sp
import math

i=0

def f(x):
    return math.sin(10*(x**0.5))**2

exact_integral= 0.5 +(1-math.cos(20))/400-math.sin(20)/20

def trapezoidal(xmin,xmax):
    return 0.5*(xmax-xmin)*(f(xmin)+f(xmax))

def adaptive_trapezoidal(xmin,xmax,whole,error):
    global i    
    left=trapezoidal(xmin,0.5*(xmin+xmax))
    right=trapezoidal(0.5*(xmin+xmax),xmax)
    delta=left+right-whole    
    if(abs(delta)<=3*error):
        i+=1
        return left+right+delta/3
    return adaptive_trapezoidal(xmin,0.5*(xmin+xmax),left,error/2)+adaptive_trapezoidal(0.5*(xmin+xmax),xmax,right,error/2)

def adaptive_trapezoid_integration(xmin,xmax,error):
    ans=adaptive_trapezoidal(xmin,xmax,trapezoidal(xmin,xmax),error)
    err=abs(ans-exact_integral)
    return ans,i,err

error=1e-6
#print(adaptive_trapezoid_integration(0,1,error))

I=[]

def trapezoid_integration(xmin,xmax,N):
    total=0    
    diff=xmax-xmin    
    for k in range(1,N):
        total+=2*f(xmin+k*diff/N)
    total+=f(xmax)+f(xmin)
    total/=(2*N)
    error=abs(total-exact_integral)
    I.append([total,N,round(error,8)])
    return total,error  

def adaptive_integration(xmin,xmax,N,error):   
    cur,delta=trapezoid_integration(xmin,xmax,N)
    if(delta<=error):
        return cur
    else:
        return adaptive_integration(xmin,xmax,2*N,error)

def adaptive_integrate(xmin,xmax,error):
    return adaptive_integration(xmin,xmax,1,error)

adaptive_integral=adaptive_integrate(0,1,error)

romberg=[]
romberg_error=[]
l=len(I)
for m in range(l):
    romberg.append([I[m][0]])
    for k in range(1,m+1):
        c=1<<(2*k)
        integral=(c*romberg[m][k-1]-romberg[m-1][k-1])/(c-1)
        romberg[m].append(integral)
    err=abs(romberg[m][m]-exact_integral)
    romberg_error.append(round(err,8))
    if(err<error):
        break

print("Adaptive Trapezoidal method [integral,n,error]\n")

for it in I:
    it[0]=round(it[0],7)
    print(it)

print("\nI = %.6f \n" % adaptive_integral)

l=len(romberg)

print("\nRomberg Table [integral] [n,error]\n")
    
for m in range(l):
    for k in range(m+1):
        romberg[m][k]=round(romberg[m][k],8)
    print(romberg[m],[1<<m,romberg_error[m]])

romberg_integral=romberg[l-1][l-1]
print("\nI = %.6f" % romberg_integral)
