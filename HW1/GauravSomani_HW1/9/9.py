import math

def primes(n):
    primes_list=[2]
    i=3
    for i in range(3,n+1):
        root_i=math.sqrt(i)
        for j in primes_list:
            if(i%j==0):
                break
            if(j>root_i):
                primes_list.append(i)
                break
    return primes_list

print(primes(10000))
