import random
from math import acos
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

theta_unicode='\u03B8'
phi_unicode='\u03D5'

def wulff_net():
    th_d = [-90,90]
    theta = np.linspace(th_d[0],th_d[1],37,dtype=int)

    def great_circle(theta):
        return -1/np.tan(theta),1/abs(np.sin(theta))
    def small_circle(theta):
        return 1/np.sin(theta),1/abs(np.tan(theta))
    
    fig, ax = plt.subplots()
    for th in theta:
    	th_r = np.deg2rad(th)
    	gc,gr = great_circle(th_r)
    	sc,sr = small_circle(th_r)
    	st = [90+th,90-th]
    	gt = [180+th,180-th]
    	if(th>0):
    		st = [-st[0],-st[1]]
    		gt = [-th,th]
    	lw=1
    	if(th%15==0):
    		lw=2
    		s = str(abs(th))+'$^o$'
    		if(th>0):
    			v = 'bottom'
    			s += 'N'
    		elif(th<0):
    			v = 'top'
    			s += 'S'
    		else:
    			v = 'center'
    		if(abs(th)!=90):
    			plt.text(np.cos(th_r), np.sin(th_r),s,horizontalalignment='left',verticalalignment=v)
    			plt.text(-np.cos(th_r), np.sin(th_r),s,horizontalalignment='right',verticalalignment=v)
    	arc1=patches.Arc((gc, 0), 2*gr, 2*gr,theta1=gt[0],theta2=gt[1],lw=lw)
    	arc2=patches.Arc((0, sc), 2*sr, 2*sr,theta1=st[0],theta2=st[1],lw=lw)
    	ax.add_artist(arc1)
    	ax.add_artist(arc2)
    plt.text(0, 1,'$90^o N$',horizontalalignment='center',verticalalignment='bottom')
    plt.text(0, -1,'$90^o S$',horizontalalignment='center',verticalalignment='top')
    plt.plot([0,0], [-1, 1], [-1, 1],[0,0], color ='black',lw=2)
    plt.title('Wulff net')
    ax.set_aspect('equal')
    plt.grid(True,'both')
    plt.xlim([-1.125,1.125])
    plt.ylim([-1.125,1.125])
    return fig,ax

def spherical_cartesian(theta,phi):
    r,z = np.sin(theta),np.cos(theta)
    x,y = r*np.cos(phi),r*np.sin(phi)
    return x,y,z

def stereo_spherical(theta,phi):
    x,y,z = spherical_cartesian(theta,phi)
    return x/(1+z),y/(1+z)

def stereo_cartesian(x,y,z):
    m = np.linalg.norm([x,y,z]) 
    return x/(m+z),y/(m+z)

def random_gen():
	phi = 2*np.pi*random.random()
	theta = acos(1-2*random.random())
	return theta,phi

#fig,ax = wulff_net()
#plt.show()
r = random_gen()
print("Randomly generated point : \n"+ theta_unicode+" = %f " %r[0])
print(phi_unicode + " = %f " %r[1])
