import numpy as np
from numpy import linalg

A = np.matrix('1 4 8 4; 4 2 3 7;8 3 6 9;4 7 9 2') 

def QR(A):
    m,n=A.shape
    Q=np.copy(A).astype(float)    
    R=np.identity(n)
    for i in range(n):
        R[i,i]=np.dot(Q[:,i],Q[:,i])**0.5
        Q[:,i] = Q[:,i]/R[i,i]
        for j in range(i+1,n):
            R[i,j] = np.dot(Q[:,i],Q[:,j])
            Q[:,j] = Q[:,j]-R[i,j]*Q[:,i]
    return Q,R

def eigen_QR(A,eps,Qk=np.identity(len(A))):
    n=A.shape[1]
    Q,R = QR(A)
    B = np.dot(R,Q)
    Qk = np.dot(Qk,Q)
    for i in range(n-1):
        if(max(abs(B[i,i+1:n]))>eps):
            return eigen_QR(B,eps,Qk=Qk)
    return np.diagonal(B),Qk.T

print("Eigenvalues and eigenvectors {Rows are eigenvectors} (by QR algorithm):")
val,vec = eigen_QR(A,1e-6)
print(val,vec)
print("\nEigenvalues (with built-in numpy eigenvals function):")
print(linalg.eigvals(A))
