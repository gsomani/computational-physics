def quadratic_1(a,b,c):
    root_D=(b*b-4*a*c)**0.5
    return (-b+root_D)/(2*a),(-b-root_D)/(2*a)

def quadratic_2(a,b,c):
    root_D=(b*b-4*a*c)**0.5
    return (2*c)/(-b-root_D),(2*c)/(-b+root_D)

def quadratic(a,b,c):
    if(b>=0):
        return quadratic_2(a,b,c)[0],quadratic_1(a,b,c)[1]
    else:
        return quadratic_1(a,b,c)[0],quadratic_2(a,b,c)[1]

a,b,c=[0.001,1000,0.001]
q1=quadratic_1(a,b,c)
q2=quadratic_2(a,b,c)
q=quadratic(a,b,c)
print("Roots using 1st method are " + str(q1))
print("Roots using 2nd method are " + str(q2))
print("Roots using composite method are " + str(q))
