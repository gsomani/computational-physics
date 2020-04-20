import math
    
approx= '\u2248'

def coeff(E,V):
    k1=math.sqrt(E)
    k2=math.sqrt(E-V)
    T=4*k1*k2/((k1+k2)*(k1+k2))
    R=1-T
    return T,R

T,R=coeff(10,9)
print("T "+ approx +" %.4f" % T)
print("R "+ approx +" %.4f" % R)
