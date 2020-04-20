import matplotlib.pyplot as plt

file = open('velocities.txt', 'r')
time=[]
vel=[]
x=[0]
N=100
for line in file: 
    t,v = line.split()
    time.append(int(t))
    vel.append(float(v))    

for i in range(N): 
    x.append(x[i]+(vel[i]+vel[i+1])/2)

plt.ylabel('Velocity(m/s)',fontsize=16)
plt.xlabel('Time(s)',fontsize=16)   
plt.grid(True, which="both")
plt.plot(time,x,label="distance-time")
plt.plot(time,vel,label="velocity-time")
plt.gca().legend()
plt.savefig('distance-time.png')
plt.show()
