import numpy as np
import matplotlib.pyplot as plt

def tdma_reduced(a,b,c):
    N=len(b)
    for i in range(N):
        b[i] -= c[i-1]*a[i]
        c[i] /= b[i]
    return b,c

def solve_tdma_reduced(a,b,c,d):
    N=len(d)
    for i in range(N):
        d[i]=(d[i]-d[i-1]*a[i])/b[i]
    for i in range(N-2,-1,-1):
        d[i] -= c[i]*d[i+1]
    return d

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def init_f(N):
    f=np.zeros(N)
    r=0.1
    xr=[round((0.5-r)*N),round((r+0.5)*N)]
    f[xr[0]+1:xr[1]] = np.ones(xr[1]-xr[0]-1)
    f[xr[0]] = xr[0] + 0.5 - (0.5-r)*N
    f[xr[1]] = 0.5 - xr[1] + (0.5+r)*N
    return f

def init(h,N):
    f=init_f(N)
    m = h*N*N
    a,b = -m*np.ones(N),(1+2*m)*np.ones(N)
    c,u = 1*a,np.zeros(N) 
    a[0],c[N-1] = 0,0    
    u[0],u[N-1] = 1,m
    v = -u
    b[0] -= u[0]*v[0]
    b[N-1] -= u[N-1]*v[N-1]    
    b_,c_=tdma_reduced(a,b,c)
    q = solve_tdma_reduced(a,b_,c_,u)
    w = 1+np.dot(v,q)
    return v,a,b_,c_,q,w,f

def fft_solution(N,t):
    f=init_f(N)
    F = np.fft.rfft(f)
    solution = np.fft.irfft(F * np.exp(-(2*np.pi*np.arange(1+N//2))**2*t))
    return solution 

def update(v,a,b,c,q,w,f):
    y = solve_tdma_reduced(a,b,c,f)
    f = y - (np.dot(v,y)/w)*q
    return f

def l1_error(f,f_true):
    M = len(f_true)
    N = len(f)
    f_ref = f_true[::M//N]
    return sum(abs(f-f_ref))/N

h=2**-14
N=2**8
t=0
f = np.empty(N+1)
v,a,b,c,q,w,f[:N] = init(h,N) 
x = np.linspace(0,1,N+1)

while(t<0.025):
    f[:N] = update(v,a,b,c,q,w,f[:N])   
    t+=h
f[N] = f[0]
plt.plot(x,f)
plt.xlabel('x')
plt.ylabel('f')
plt.grid(True,'both')
plt.title('f(t=0.025) {256 points}')
plt.savefig('diffusion_256.png')
plt.show()

N=2**7
l1,dx=np.empty(3),np.empty(3)
h=2**-14
g = []
for i in range(4):
	t=0
	f = np.empty(N+1)
	v,a,b,c,q,w,f[:N] = init(h,N)
	while(t<0.25):
		f[:N] = update(v,a,b,c,q,w,f[:N])
		t+=h
	x = np.linspace(0,1,N+1)
	f[N]=f[0]
	plt.plot(x,f,label=str(N)+'points')
	g.append(f[:N])
	N*=2

plt.xlabel('x')
plt.ylabel('f')
plt.grid(True,'both')
plt.gca().legend()	
plt.title('Convergences test for f(t=0.25) ')
plt.savefig('diffusion_convergence.png')
plt.show()

for i in range(3):
    l1[i] = l1_error(g[i],g[3])
    dx[i] = 1/len(g[i])

slope = 0.5*np.log2(l1[0]/l1[2])
print("Slope of L1-error plot is %.1f." % slope)
order = round(slope)
print("Order of convergence is %i." % order)
plt.xlabel('$\Delta x$')
plt.ylabel('L1 error')
plt.loglog(dx,l1,basex=2,basey=2)
plt.grid(True,'both')
plt.title('L1 Error')
plt.savefig('error.png')
plt.show()
