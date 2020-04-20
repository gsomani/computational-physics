def float_factorial(n):
    if(n==0):
        return 1
    return float(n)*float_factorial(n-1)

def int_factorial(n):
    if(n==0):
        return 1
    return n*int_factorial(n-1)

print("Factorial (using integer) = %i" %int_factorial(200))
print("Factorial (using float) = %f" %float_factorial(200))
