import scipy as sp
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,15)

file = open('stm.txt', 'r')
height=sp.empty([663, 676])
i=0
for line in file:
    height[i]=line.split()
    i+=1
plt.imshow(height)
plt.title('(111) Silicon surface mapping using STM') 
plt.colorbar().set_label('Height (STM data)')
plt.savefig('silicon.png')
plt.show()
