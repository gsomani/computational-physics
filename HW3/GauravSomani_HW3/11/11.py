def f(x,y):
    return y*y*(1-x)-x*x*x

def g(x,y):
    return y*y + x*x -1

def fx(x,y):
    return -3*x*x-y*y

def fy(x,y):
    return 2*y*(1-x)

def gx(x,y):
    return 2*x

def gy(x,y):
    return 2*y

def newton_iteration(x,y):
    det=fx(x,y)*gy(x,y)-gx(x,y)*fy(x,y)
    diff_x=(fy(x,y)*g(x,y)-gy(x,y)*f(x,y))/det
    diff_y=(gx(x,y)*f(x,y)-fx(x,y)*g(x,y))/det
    return x+diff_x,y+diff_y

def newton(guess,eps):
    x,y=newton_iteration(guess[0],guess[1])
    diff=[abs(x-guess[0]),abs(y-guess[1])]
    diff_size=diff[0]*diff[0]+diff[1]*diff[1]
    if(diff_size<0.25*eps*eps):
        return x,y
    return newton([x,y],eps)

root=[[0.6,-0.8],[0.6,0.8],[-1.6,1.25j],[-1.6,-1.25j]] 
eps=1e-6

print("Roots : \n")
for i in range(2):
    root[i]=newton(root[i],eps)
    print("x = %.6f, y = %.6f" %(root[i][0],root[i][1]))
for i in range(2,4):
    root[i]=newton(root[i],eps)
    print("x = %.6f, y = %.6f i" %(root[i][0].real,root[i][1].imag))
