import scipy as sp
from scipy import constants
import math
import time
    
approx= '\u2248'
G=constants.G
pi=constants.pi
sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

T_day=3600*24
M_sun=1.989e30

def orbit(l1,v1):
    b=2*G*M_sun/(v1*l1)
    c=v1*v1-2*G*M_sun/l1
    D=b*b+4*c
    v2=(b-math.sqrt(D))/2
    l2=(v1*l1)/v2
    e=(l2-l1)/(l1+l2)
    T=pi*(l2+l1)*(l2+l1)*math.sqrt(1-e*e)/(2*v1*l1)
    T=T/T_day # to get T in days
    return l2,v2,T,e

l1=float(input("l1".translate(sub)+" (in m) = "))
v1=float(input("v1".translate(sub)+" (in m/s) = "))
l2,v2,T,e=orbit(l1,v1)

print("l2 ".translate(sub)+ approx +" %.4E m" % l2)
print("v2 ".translate(sub)+ approx +" %.4E m/s" % v2)
print("T "+ approx +" %.2f days" % T)
print("e "+ approx +" %.4f \n" % e)

l1_earth=1.471e11
v1_earth=3.0287e4
l1_comet=8.783e10
v1_comet=5.4529e4
l2_earth,v2_earth,T_earth,e_earth=orbit(l1_earth,v1_earth)
l2_comet,v2_comet,T_comet,e_comet=orbit(l1_comet,v1_comet)

print("Orbit parameters for Earth :")
print("l2 ".translate(sub)+ approx +" %.4E m" % l2_earth)
print("v2 ".translate(sub)+ approx +" %.4E m/s" % v2_earth)
print("T "+ approx +" %.2f days" % T_earth)
print("e "+ approx +" %.4f \n" % e_earth)

print("Orbit parameters for Halley's comet :")
print("l2 ".translate(sub)+ approx +" %.4E m" % l2_comet)
print("v2 ".translate(sub)+ approx +" %.4E m/s" % v2_comet)
print("T "+ approx +" %.2f days" % T_comet)
print("e "+ approx +" %.4f \n" % e_comet)
