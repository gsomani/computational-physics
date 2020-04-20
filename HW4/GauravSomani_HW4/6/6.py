import numpy as np
from cmath import exp
from math import pi

def fft(y):
    n=len(y)
    if(n==1):
        return y    
    y0=y[0:n:2]
    y1=y[1:n:2]
    c0=fft(y0)
    c1=fft(y1)
    num=-2j*pi/n
    k=n//2
    c=np.empty(n,dtype=complex)
    for i in range(k):
        a=exp(i*num)
        c[i]=c0[i]+a*c1[i]
        c[i+k]=c0[i]-a*c1[i]
    return c

file = open('pitch.txt', 'r')
pitch=[]

for line in file:
    p=float(file.readline())
    pitch.append(p)

np_fft = np.fft.fft(np.array(pitch))
cooley_fft = fft(pitch)

print("FFT (calculated by NumPy)")
print(np_fft)
print("FFT (calculated by function defined in this program)")
print(cooley_fft)
print("Magnitude of difference (at each frequency) between both methods")
print(abs(np_fft-cooley_fft))
