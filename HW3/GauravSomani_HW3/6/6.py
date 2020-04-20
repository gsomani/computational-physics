import scipy.constants as constants
import math

h=constants.h
c=constants.c
k=constants.k
b1=h*c/k

def f(x):
    return 5*math.exp(-x)+x-5

left_limit=4.5
right_limit=5

def bisection(left,right,error):
    mid=(left+right)/2
    if(abs(mid-left)<error):
        return mid
    if(f(mid)*f(left)<0):
        return bisection(left,mid,error)
    return bisection(mid,right,error)

error=1e-6
x=bisection(left_limit,right_limit,error)
b=b1/x
print("b = %.6f m-K" %b)
lamda=5.02e-7
T=b/lamda
print("T = %.1f K" %T)



