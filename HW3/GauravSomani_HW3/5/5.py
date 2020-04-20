def relaxation(a,b,x,eps):
    x0=x    
    y=b/(a+x*x)
    x=y*(a+x*x)
    diff=x-x0
    if(abs(diff)<eps):
        return x,y
    return relaxation(a,b,x,eps)

print("Relaxation  method:\n")
print("(x,y) = "+ str(relaxation(1,2,1,1e-6)))
