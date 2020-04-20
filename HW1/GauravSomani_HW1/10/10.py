import time

sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

def catalan_recursive(n):
    if(n==0):
        return 1
    return catalan_recursive(n-1)*(4*n-2)//(n+1)

start_time = time.time()
print("C100".translate(sub)+" = %i" % catalan_recursive(100))
print("Time taken = %.6f s" % (time.time()-start_time))

def gcd(m,n):
    if(n==0):
        return m
    return int(gcd(n,m%n))
g=gcd(108,192)
print("g(108,192) = %i" % gcd(108,192))
