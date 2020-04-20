import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def load_music(filename):
    M=10000
    rate=44100
    file = open(filename, 'r')
    music=[]

    for line in file:
        p=float(file.readline())
        music.append(p)

    plt.grid(True, which="both")
    plt.xlabel('n')
    plt.ylabel('$y_n$')
    instrument=filename[:len(filename)-4]
    plt.title('Waveform ('+instrument+')')     
    plt.plot(music)
    plt.savefig('waveform ('+instrument+')'+'.png')    
    plt.show()

    c = np.fft.rfft(music)
    A = abs(c)
    N = len(music)

    plt.grid(True, which="both")    
    plt.xlabel('k')
    plt.ylabel('$|c_k|$')
    plt.title('Discrete Fourier transform ('+instrument+')')     
    plt.plot(A[0:M])
    plt.savefig('fft ('+instrument+')'+'.png')    
    plt.show()
    
    k=np.argmax(A[0:M])
    return k*rate/N

f_piano=load_music('piano.txt')
f_trumpet=load_music('trumpet.txt')
f_middleC=261
ratio_piano=round(math.log2(f_piano/f_middleC))
ratio_trumpet=round(math.log2(f_trumpet/f_middleC))
print("Note played by piano is %i octaves above middle C." %ratio_piano)
print("Note played by trumpet is %i octaves above middle C." %ratio_trumpet) 
