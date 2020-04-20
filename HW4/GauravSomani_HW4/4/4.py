import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

file = open('dow.txt', 'r')
price=[]
for line in file:
    p=float(line)
    price.append(p)

plt.grid(True, which="both")
plt.plot(price)
plt.xlabel('Day')
plt.ylabel('Price')
plt.xticks(np.linspace(0,1023,17))
plt.title('Dow Jones Industrial Average (Raw data)') 
plt.savefig('dow_data.png')       
plt.show()

c = np.fft.rfft(price)
N=len(c)

for i in range(N//10,N):
    c[i]=0

y = np.fft.irfft(c)

plt.grid(True, which="both")
plt.plot(price,color='blue',label='Original data')
plt.plot(y,color='red',label='10% fourier coefficients with FFT')
plt.title('Dow Jones Industrial Average (10% fourier coefficients with FFT)')     
plt.gca().legend()
plt.savefig('dow_fft_10percent.png')
plt.show()

for i in range(N//50,N//10):
    c[i]=0

y = np.fft.irfft(c)

plt.grid(True, which="both")
plt.plot(price,color='blue',label='Original data')
plt.plot(y,color='red',label='2% fourier coefficients with FFT')
plt.title('Dow Jones Industrial Average (2% fourier coefficients with FFT)')
plt.gca().legend()
plt.savefig('dow_fft_2percent.png')    
plt.show()

M=1000

ys=np.empty(M)

for i in range(M):
    ys[i]= -1+2*(i<(M//2))    

cs = np.fft.rfft(ys)

for i in range(10,1+M//2):
    cs[i]=0

yr = np.fft.irfft(cs)

plt.grid(True, which="both")
plt.plot(ys,color='blue',label='Original data')
plt.plot(yr,color='red',label='10 fourier coefficients with FFT')
plt.title('Square Wave (10 fourier coefficients with FFT)')
plt.yticks(np.linspace(-1.2,1.2,13))
plt.gca().legend()
plt.savefig('square_wave_fft.png')
plt.show()
