import scipy as sp
import matplotlib.pyplot as plt

file = open('sunspots.txt', 'r')
month=[]
sunspot=[]
for line in file: 
    m,s = line.split()
    month.append(int(m))
    sunspot.append(float(s))    

plt.ylabel('Sunspots',fontsize=16)
plt.xlabel('Month',fontsize=16)        
plt.xticks(sp.linspace(0, 3250, 14))
plt.yticks(sp.linspace(0, 260,27))
plt.grid(True, which="both")
plt.plot(month,sunspot)
plt.title('(a)')
plt.savefig('all_months.png')    
plt.show()

plt.ylabel('Sunspots',fontsize=16)
plt.xlabel('Month',fontsize=16)        
plt.xticks(sp.linspace(0, 1000, 11))
plt.yticks(sp.linspace(0, 260, 27))
plt.grid(True, which="both")
plt.plot(month[0:1000],sunspot[0:1000])
plt.title('(b)')
plt.savefig('1000_months.png')    
plt.show()
         
sunspot_avg=sp.array(sunspot)
months=len(month)

for i in range(6):
    sunspot_avg[0]+=sunspot[i]
    
for i in range(1,1000):
    subtract=0
    if(i-6>=0):
        subtract=sunspot[i-6]    
    sunspot_avg[i]=sunspot_avg[i-1]+sunspot[i+5]-subtract

for i in range(1000):
    sunspot_avg[i]/=(6+min(i,5))

plt.ylabel('Sunspots',fontsize=16)
plt.xlabel('Month',fontsize=16)        
plt.xticks(sp.linspace(0, 1000, 11))
plt.yticks(sp.linspace(0, 260, 27))
plt.grid(True, which="both")
plt.title('(c)')
plt.plot(month[0:1000],sunspot[0:1000],label='Original data')
plt.plot(month[0:1000],sunspot_avg[0:1000],label='Running average')
plt.gca().legend()
plt.savefig('running_average.png')
plt.show()
