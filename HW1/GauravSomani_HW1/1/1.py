import scipy as sp
from scipy import constants
import math
    
approx= '\u2248'
G=constants.G
pi=constants.pi

M_earth=5.974e24
R_earth=6.378e6

def altitude(T):
    r=(G*M_earth*T*T/(4*pi*pi))**(1/3)
    return r-R_earth

def print_altitude(T):
    H = altitude(T)
    H_km = H/1000
    if(H<0):
        print("Not possible\n")
    else:
        print("Altitude for orbit period of %i seconds " %T + approx +" %.4E m " %H + approx +" %.2f km\n" %H_km)

print("(a)")        
T = float(input('T (in seconds) = '))
print_altitude(T)
        
T_minute=60
T_hour=3600
T_day=T_hour*24

print("(b)")        
print("Satellites oribiting once a day :")
print_altitude(T_day)      
print("Satellites oribiting once every 90 minutes :")
print_altitude(T_minute*90)
print("Satellites oribiting once every 45 minutes :")
print_altitude(T_minute*45)

print("(c)")        
print("24 hours is one solar day while 23.93 hours is one sidereal day. Solar day means time period between two consecutive events of facing the sun for a point on the earth. Sidereal day is the period of rotation of earth about its axis passing through center of earth. Due to simulaneous rotation and revolution of earth, after one sidereal day, a point on the earth which was facing the sun at the start of the day does not face sun. Earth needs to rotate for some more time for that point to face the sun. So, solar day is longer than earth's rotation period about its axis. But geostationary satellite needs to be pointing to a particular point on earth. So, it should rotate with that point. Hence, it should orbit once in sidereal day.\n")
T_geosynchronous=T_hour*23.93
print("Geosynchronous satellites:")
print_altitude(T_geosynchronous)
diff_km=(altitude(T_day)-altitude(T_geosynchronous))/1000
print("Difference in altitude is about %.2f km. Geostationery satellite needs to be about %.2f km lower than satellite with orbit period 24 hours." % (diff_km,diff_km))   

