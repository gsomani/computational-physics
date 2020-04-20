import math
from scipy import constants

approx= '\u2248'
theta_unicode='\u03B8'
pi=constants.pi

def cartesian_to_polar(x,y):
    r=math.sqrt(x*x+y*y)
    if(x!=0):
        theta=180*math.atan(y/x)/pi
        if(x<0):
            theta-=180
        if(theta<-180):
            theta+=360
    else:
        if(y>0):
            theta=90
        else: 
            theta=-90
    return r,theta

x = float(input('x = '))
y = float(input('y = '))
r,theta=cartesian_to_polar(x,y)
print("r,"+theta_unicode +' ' + approx + " ( %.4f , %.4f )" % (r,theta))
