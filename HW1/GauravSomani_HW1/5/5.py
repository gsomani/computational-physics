import math
import time

sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    
approx= '\u2248'

def catalan(n):
    c=1
    for i in range(n):
        c=c*(4*i+2)//(i+2)
    return c    

start_time = time.time()
print("C100".translate(sub)+" = %i" % catalan(100))
print("Time taken "+ approx +" %.6f s" % (time.time()-start_time))
