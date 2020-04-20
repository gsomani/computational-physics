import scipy as sp

n=200
root=n**0.5
B=2*sp.diag(sp.full(n,1.0))+sp.random.rand(n,n)/(n**0.5)
A=0.5*(B+sp.transpose(B))
b=sp.full(n,1.0)    

def solve_triangular(A,b,c):
    y=sp.copy(b)
    m,n=A.shape    
    if(c=='l'):
        for i in range(m):
            y[i] = ( y[i] - sp.dot(A[i,0:i],y[0:i]) ) / A[i,i]
    else:
        for i in range(n-1,-1,-1):
            y[i] = ( y[i] - sp.dot(A[i,i+1:n],y[i+1:n]) ) / A[i,i]
    return y[0:n]
    
def solve_lu(A,b):
    m,n=A.shape
    U=sp.copy(A).astype(float)
    L=sp.identity(m)
    for k in range(m-1):
        i=k+sp.argmax(abs(U[k:,k]))
        s=U[[i,k]]
        U[k,k:m],U[i,k:m]=s[0,k:m],s[1,k:m]
        s=L[[i,k]]
        L[k,0:k],L[i,0:k]=s[0,0:k],s[1,0:k]        
        b[[k,i]]=b[[i,k]]
        for j in range(k+1,m):
            L[j,k]=U[j,k]/U[k,k]
            U[j,k:m]=U[j,k:m]-L[j,k]*U[k,k:m]
    y=solve_triangular(L,b,'l')    
    x=solve_triangular(U,y,'u')
    return L,U,x

def QR(A):
    m,n=A.shape
    Q=sp.copy(A).astype(float)
    R=sp.identity(n)
    for i in range(n):
        R[i,i]=sp.dot(Q[:,i],Q[:,i])**0.5
        Q[:,i] = Q[:,i]/R[i,i]
        for j in range(i+1,n):
            R[i,j] = sp.dot(Q[:,i],Q[:,j])
            Q[:,j] = Q[:,j]-R[i,j]*Q[:,i]
    return Q,R

def solve_QR(A,b):
    Q,R = QR(A)
    y=sp.dot(sp.transpose(Q),b)
    x=solve_triangular(R,y,'u')
    return x

L,U,x_lu=solve_lu(A,b)
x_qr=solve_QR(A,b)

print("x (by LU decomposition)")
print(x_lu)
print("\nx (by QR decomposition)")
print(x_qr)
