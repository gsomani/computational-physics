import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

N=1000
pi=np.pi

y = [np.empty(N), np.arange(N), np.empty(N)]

for i in range(N):
    y[0][i] = -1+2*(i<(N//2))    
    y[2][i] = math.sin(pi*i/N)*math.sin(20*pi*i/N)

title = ['Square Wave','Sawtooth wave', 'Modulated sine wave']
filename = ['square.png','sawtooth.png','modulated_sine.png']

for i in range(3):
    c = np.fft.rfft(y[i])
    plt.grid(True, which="both")
    plt.plot(abs(c))
    plt.title(title[i])
    plt.xlabel('k')
    plt.ylabel('$|c_k|$') 
    plt.savefig(filename[i])
    plt.show() 
