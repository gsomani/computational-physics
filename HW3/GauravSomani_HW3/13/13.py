import scipy as sp
from scipy import linalg

A = sp.matrix('1 4 8 4; 4 2 3 7;8 3 6 9;4 7 9 2') 

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

def eigen_QR(A,eps):
    n=A.shape[1]
    Q,R = QR(A)
    B = sp.dot(R,Q)
    for i in range(n-1):
        if(max(abs(B[i,i+1:n]))>eps):
            return eigen_QR(B,eps)
    return sp.diagonal(B)

print("Eigenvalues (by QR algorithm):")
print(eigen_QR(A,1e-6))
print("\nEigenvalues (with built-in scipy eigenvals function):")
print(linalg.eigvals(A).real)
