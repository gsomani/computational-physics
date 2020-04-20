approx= '\u2248'

def binomial(n,k):
    if(k==0):
        return 1
    return int(binomial(n,k-1)*(n-k+1)/k)

def pascal_triangle_line(n):
    coeff_list=[]
    for k in range(n+1):
        coeff_list.append(binomial(n,k))
    return coeff_list

def print_pascal_triangle(n):
    for i in range(n):
        print(pascal_triangle_line(i+1))

print_pascal_triangle(20)

n=100
heads=60
P_60=binomial(100,60)/(1<<n);
P_more_than_60=0
for k in range(heads,n+1):
    P_more_than_60+=binomial(n,k)
P_more_than_60/=1<<n

print("Probability of 60 heads "+ approx +" %.4f" % P_60)
print("Probability of more than 60 heads "+ approx +" %.4f" % P_more_than_60)
