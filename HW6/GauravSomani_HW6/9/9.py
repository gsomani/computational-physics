import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

N=40
L=2
r=0.3

def rho(y,x):
    if(abs(x)<r):
        return np.exp(-abs(y)/0.1)
    return 0
    
def jacobi(f,N,L):
    cur=np.zeros([N+1,N+1])
    for i in range(1,N):
        for j in range(1,N):
            cur[i][j] = (f[i-1][j]+f[i+1][j]+f[i][j-1]+f[i][j+1] - rho(-1+i*L/N,-1+j*L/N)*L*L/(N*N))/4    
    return cur,np.sum(abs(cur-f))/((N-1)*(N-1))

def sor(f,N,L,w):
    g = np.copy(f)    
    for i in range(1,N):
        for j in range(1,N):    
            r = f[i-1][j]+f[i+1][j]+f[i][j-1]+f[i][j+1] - 4*f[i][j] - rho(-1+i*L/N,-1+j*L/N)*L*L/(N*N)
            f[i][j] += w*r/4
    return f,np.sum(abs(g-f))/((N-1)*(N-1))

it=1000
phi_j = np.zeros([N+1,N+1])
err = np.zeros([5,it])
x = np.linspace(-L/2,L/2,N+1)

for i in range(1,it+1):
    phi_j,err[0][i-1]=jacobi(phi_j,N,L)
    if(i==10 or i==100 or i==1000):
        plt.contourf(x,x,phi_j)
        plt.gca().set_aspect('equal')
        plt.colorbar().set_label('Potential')
        plt.title('Potential {Jacobi}(iterations = ' + str(i)+' )')
        plt.savefig('potential_jacobi'+str(i)+'.png')    
        plt.show()

plt.loglog(err[0])
plt.xlim([1,it])
plt.grid(True,'both') 
plt.ylabel('L1 Error')
plt.xlabel('Number of iterations')
plt.title('L1 error (Jacobi method)')
plt.savefig('error_jacobi.png') 
plt.show()

w_opt = 2/(1+np.pi/N)
w = [1,0.5,1.5,w_opt]

phi = np.zeros([4,N+1,N+1])

for j in range(1,5):
    for i in range(1,it+1):    
        phi[j-1],err[j][i-1]=sor(phi[j-1],N,L,w[j-1])
        if(i==10 or i==100 or i==1000):
            plt.contourf(x,x,phi[j-1])
            plt.gca().set_aspect('equal')
            plt.colorbar().set_label('Potential')
            if(w[j-1]==1):
                plt.title('Potential {Gauss-Seidel} (iterations = ' + str(i)+' )')
                plt.savefig('potential_gs'+str(i)+'.png')
            else:    
                plt.title('Potential {SOR, $\omega$ = '+str(round(w[j-1],2))+' } (iterations = ' + str(i)+' )')
                plt.savefig('potential_sor'+str(i)+'.png')
            plt.close()            
    
for j in range(1,5):
    if(w[j-1]==1):
        plt.loglog(err[j],label='Gauss-Seidel')
    else:    
        plt.loglog(err[j],label='SOR ( $\omega$ = '+str(round(w[j-1],2))+' )')

plt.gca().legend()
plt.grid(True,'both') 
plt.xlim([1,it])
plt.ylabel('L1 Error')
plt.xlabel('Number of iterations')
plt.title('L1 error (SOR and Gauss-Siedel)')
plt.savefig('error_sor.png') 
plt.show()

plt.loglog(err[0],label='Jacobi')
plt.loglog(err[1],label='Gauss-Seidel')
for j in range(2,5):
    plt.loglog(err[j],label='SOR ( $\omega$ = '+str(round(w[j-1],2))+' )')
plt.gca().legend()
plt.grid(True,'both')
plt.xlim([1,it]) 
plt.ylabel('L1 Error')
plt.xlabel('Number of iterations')
plt.title('L1 error (Comparision of convergence rates of different methods)')
plt.savefig('error_comparision.png') 
plt.show()
