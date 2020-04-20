import numpy as np
import matplotlib.pyplot as plt

def p(t,tau):
    return 1-2**(-t/tau)    

def Pb_decay(x,p):
    r = np.sum( np.random.rand(x[1]) <= p)
    x[1] -= r
    x[0] += r
    return x

def Tl_decay(x,p):
    r = np.sum( np.random.rand(x[2]) <= p)
    x[2] -= r
    x[1] += r
    return x

def Bi_decay(x,p,route_prob):
    r = np.sum( np.random.rand(x[3]) <= p)
    r_dec = np.sum( np.random.rand(r) <= route_prob)   
    return x[0],x[1]+r-r_dec,x[2]+r_dec,x[3]-r

def decay(X,dt,time,tau,Tl_prob):
    N = round(time/dt)
    x = np.empty([N+1,4],dtype=int)
    x[0] = X
    prob = np.empty(3)
    for i in range(3):
         prob[i] = p(dt,tau[i])
    for t in range(N):
       cur = Pb_decay(x[t],prob[0])
       cur = Tl_decay(cur,prob[1])
       x[t+1] = Bi_decay(cur,prob[2],Tl_prob) 
    return x,np.linspace(0,time,N+1)

time = 20000
dt = 1
X = [0,0,0,10000]
tau = [198,132,2760]
Tl_prob = 0.0209
x,t = decay(X,dt,time,tau,Tl_prob)
isotope = ["$^{209}$Bi","$^{209}$Pb","$^{209}$Ti","$^{213}$Bi"]
for i in range(4):
    plt.plot(t,x[:,i],label=isotope[i])
plt.gca().legend()
plt.title('Radioactive decay chain')
plt.xlabel('time (s)')
plt.ylabel('Number of atoms')
plt.grid(True, which="both")
plt.savefig('decay.png')  
plt.show()
