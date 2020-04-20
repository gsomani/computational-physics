import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (12,12)

file = open('sunspots.txt', 'r')
month=[]
sunspot=[]
for line in file: 
    m,s = line.split()
    month.append(int(m))
    sunspot.append(float(s))    

plt.ylabel('Sunspots',fontsize=16)
plt.xlabel('Month',fontsize=16)        
plt.xticks(np.linspace(0, 3250, 14))
plt.yticks(np.linspace(0, 260,27))
plt.grid(True, which="both")
plt.plot(month,sunspot)
plt.title('Sunspots')
plt.savefig('all_months.png')    
plt.show()

c = np.fft.rfft(sunspot)
ps = abs(c)**2
N=len(month)

plt.grid(True, which="both")
plt.plot(ps)
plt.title('Power spectrum')
plt.xlabel('k')
plt.ylabel('$|c(k)|^2$') 
plt.savefig('fft.png')  
plt.show()

k=np.argmax(ps[1:N])+1
period=N/k
print("k = %i , Period = %.2f months" % (k,period))
