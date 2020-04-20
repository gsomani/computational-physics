import scipy as sp
import matplotlib.pyplot as plt
from scipy import constants

approx= '\u2248'
h=constants.h
q=1.602e-19

file = open('millikan.txt', 'r')
freq=[]
voltage=[]
Exy=Exx=Ex=Ey=0
for line in file: 
    x,y = line.split()
    x=float(x)
    y=float(y)
    freq.append(x)
    voltage.append(y) 
    Exx+=x*x
    Exy+=x*y
    Ex+=x
    Ey+=y

n=len(freq)
Exx/=n
Exy/=n
Ex/=n
Ey/=n
m=(Exy-Ex*Ey)/(Exx-Ex*Ex)
c=(Exx*Ey-Ex*Exy)/(Exx-Ex*Ex)
print("m "+approx+" %.4E Vs" % m)
print("c "+approx+" %.4f V" % c) 
fitted_voltage=[]
for i in range(n):
    fitted_voltage.append(m*freq[i]+c)
plt.ylabel('Voltage (V)',fontsize=16)
plt.xlabel('Frequency (Hz)',fontsize=16)  
plt.grid(True, which="both")
plt.plot(freq,voltage,'ko',freq,fitted_voltage)
plt.title('Least-squares fitting for photoelectric effect data')
plt.savefig('line_fit.png')
plt.show()

h_calc=m*q
h_error_percent=100*(1-h_calc/h)
print("h (calculated) "+ approx +" %.4E Js" % h_calc)
print("h "+ approx +" %.4E Js" % h)
print("Error "+ approx +" %.3f" % h_error_percent+" %")
