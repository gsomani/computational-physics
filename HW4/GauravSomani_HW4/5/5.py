import numpy as np
import matplotlib.pyplot as plt
import math
import dcst

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

file = open('dow2.txt', 'r')
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

for i in range(N//50,N):
    c[i]=0

y = np.fft.irfft(c)

plt.grid(True, which="both")
plt.plot(price,color='blue',label='Original data')
plt.plot(y,color='red',label='2% fourier coefficients with FFT')
plt.title('Dow Jones Industrial Average (2% fourier coefficients with FFT)')     
plt.gca().legend()
plt.savefig('dow_fft_2percent.png')
plt.show()

dct = dcst.dct(price)
N=len(price)

for i in range(N//50,N):
    dct[i]=0

y = dcst.idct(dct)

plt.grid(True, which="both")
plt.plot(price,color='blue',label='Original data')
plt.plot(y,color='red',label='2% fourier coefficients with DCT')
plt.title('Dow Jones Industrial Average (2% fourier coefficients with DCT)')     
plt.gca().legend()
plt.savefig('dow_dct_2percent.png')
plt.show()
