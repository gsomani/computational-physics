import random
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})
plt.rcParams["figure.figsize"] = (16,16)

def random_dir():
    r = random.randint(0,3)
    a,b = r//2 , 2*(r%2)-1
    if(a):
        return 0,b
    return b,0

def brownian_update(x,y,L):
    diff = random_dir()
    x_,y_ = x+diff[0],y+diff[1]
    if(x_<0 or x_>=L or y_<0 or y_>=L):
    	return x,y
    return x_,y_

def brownian_motion(L,n):
	x = np.empty([n+1,2],dtype=int)
	x[0] = L//2,L//2
	cur=10
	for i in range(n):
		x[i+1] = brownian_update(x[i][0],x[i][1],L)
		if(i==cur-1):
			plt.plot(x[:cur,0],x[:cur,1])
			plt.title('Random Walk (n = %i)' %cur)
			plt.xlim([0,L-1])
			plt.ylim([0,L-1])
			plt.gca().set_aspect('equal')
			plt.savefig('random_walk_%07d.png' %cur)
			plt.close()
			if(i<100):
				cur *= 10
			else:
				cur *= 2
	plt.plot(x[:,0],x[:,1])
	plt.title('Random Walk (n = %i)' %n)
	plt.xlim([0,L-1])
	plt.ylim([0,L-1])
	plt.gca().set_aspect('equal')
	plt.savefig('random_walk_%i.png' %n)
	plt.show()
	return x

L = 101
n = 1000000
x = brownian_motion(L,n)
