import scipy as sp
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (15,15)

file = open('altitude.txt', 'r')
h=3e4
x=sp.linspace(0, 1023*h, 1024)
y=sp.linspace(0, 511*h, 512)
w=sp.empty([512,1024])
wx=sp.empty([512,1024])
wy=sp.empty([512,1024])
i=0
for line in file:
    w[i]=line.split()
    i+=1

def dwx(x,y,h,xmax):
    if(x==0):
        return (w[y][x+1]-w[y][x])/h
    if(x==xmax):
        return (w[y][x]-w[y][x-1])/h
    return 0.5*(w[y][x+1]-w[y][x-1])/h

def dwy(x,y,h,ymax):
    if(y==0):
        return (w[y+1][x]-w[y][x])/h
    if(y==ymax):
        return (w[y][x]-w[y-1][x])/h
    return 0.5*(w[y+1][x]-w[y-1][x])/h

for i in range(512):
    for j in range(1024):
        wx[i][j]=dwx(j,i,h,1023)
        wy[i][j]=dwy(j,i,h,511)

def I(x,y):
    wxx=wx[y][x]
    wyy=wy[y][x]
    return (wxx+wyy)/(2*(wxx*wxx+wyy*wyy+1))**0.5           

intensity=sp.empty([512,1024])

for i in range(512):
    for j in range(1024):
        intensity[511-i][j]=I(j,i)

plt.pcolormesh(x,y,intensity,cmap='gray')
plt.gca().set_aspect('equal')
plt.title('World Relief')
plt.savefig('world_relief.png')
plt.show()

file = open('stm.txt', 'r')
h=2.5
wx=sp.empty([663,676])
wy=sp.empty([663,676])
w=sp.empty([663, 676])
i=0
for line in file:
    w[i]=line.split()
    i+=1

for i in range(663):
    for j in range(676):
        wx[i][j]=dwx(j,i,h,675)
        wy[i][j]=dwy(j,i,h,662)

intensity=sp.empty([663,676])
x=sp.linspace(0, 675*h, 676)
y=sp.linspace(0, 662*h, 663)

for i in range(663):
    for j in range(676):
        intensity[662-i][j]=I(j,i)

plt.pcolormesh(x,y,intensity,cmap='gray')
plt.gca().set_aspect('equal')
plt.title('Silicon (111)')
plt.savefig('silicon.png')
plt.show()
