import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (20,20)

G=1
M = [150,200,250]

def f(r,t):
    N=len(r)//2
    r=np.array(r)
    R=np.zeros([N,N,2])
    a=np.zeros([N,2])
    for i in range(N):
        for j in range(i):
            R[i][j]=r[i]-r[j]
            R[j][i]=-R[i][j]
    for i in range(N):
        for j in range(i):
            d = la.norm(R[i][j])
            div = d*d*d
            a[i] += G*M[j]*R[j][i]/div
            a[j] += G*M[i]*R[i][j]/div
    ret=np.empty([2*N,2])
    ret[0:N],ret[N:2*N] = r[N:2*N],a   
    return ret

def rk4_step(x,t,h,f):
	k0=f(x,t)
	k1=f(x+0.5*h*k0,t+0.5*h)
	k2=f(x+0.5*h*k1,t+0.5*h)
	k3=f(x+h*k2,t+h)
	dx = h*(k0+2*k1+2*k2+k3)/6
	return dx

def adaptive_rk4(r0,v0,time,h,delta,f,mult):
    d=len(r0)
    x0 = np.empty([2*d,2])
    x0[0:d],x0[d:2*d] = r0,v0
    x,t=x0,0
    r,v,T=[x0[0:d]],[x0[d:2*d]],[t]
    cur_x=np.empty(d)
    cur_y=np.empty(d)
    for j in range(d):
        cur_x[j],cur_y[j]= r[0][j]
    X,Y=[cur_x.tolist()],[cur_y.tolist()]
    i=1
    eps=2**-52
    while(1):
        dxh = rk4_step(x,t,h,f)
        dx1 = dxh + rk4_step(x+dxh,t+h,h,f)
        dx2 = rk4_step(x,t,2*h,f)
        diff = dx1 - dx2
        diff_r = la.norm(diff[0:d])
        m=(h*delta/diff_r)**(1/4)
        if(m < 1):
            h *= m 
            continue
        T.append(t+h)
        t += 2*h
        T.append(t)    
        xh = x + dxh        
        x += dx1 + diff/15
        r.append(xh[0:d])	
        v.append(xh[d:2*d])
        r.append(x[0:d])	
        v.append(x[d:2*d])
        for k in range(2):
            for j in range(d):
                cur_x[j],cur_y[j] = r[i][j]
        X.append(cur_x.tolist())
        Y.append(cur_y.tolist())
        diff_t=time-t
        if(diff_t<eps):
            break
        h*=min(mult,m)
        i+=2
        if(diff_t<2*h):
            h=diff_t/2       
    return np.array(X),np.array(Y),T

r0 = [[3,1],[-1,-2],[-1,1]]
v0 = [[0,0],[0,0],[0,0]]
delta=2**-10
time=2

x,y,t=adaptive_rk4(r0,v0,time,0.1,delta,f,1.1)
for i in range(3):
    plt.plot(x[:,i],y[:,i],label="Star "+str(i+1))
plt.gca().set_aspect('equal')
plt.gca().legend()
plt.grid(True, which="both")
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.savefig('trail.png')
plt.show()
