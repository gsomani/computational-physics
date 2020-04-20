import math
    
approx= '\u2248'

M=0
L=600
root_2=math.sqrt(2)
root_3=math.sqrt(3)
for i in range(1,L+1):
    cur=(2/i)*(4/root_3+3)
    if(i%2==0):
        M+=cur
    else:
        M-=cur
    M+=12/(i*root_2)    
    for j in range(1,i):
        cur=24/math.sqrt(2*i*i+j*j)
        if(j%2==0):
            M+=cur
        else:
            M-=cur 
        cur=24/math.sqrt(i*i+2*j*j)
        if(i%2==0):
            M+=cur
        else:
            M-=cur
        cur=24/math.sqrt(i*i+j*j)
        if((i+j)%2==0):
            M+=cur
        else:
            M-=cur         
        for k in range(1,j):
            cur=48/math.sqrt(i*i+j*j+k*k)
            if((i+j+k)%2==0):
                M+=cur
            else:
                M-=cur 
print("M "+approx+" %.3f" % M)
