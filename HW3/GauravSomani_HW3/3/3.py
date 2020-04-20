def f(x):
    return x*(x-1)

def deriv_f(x):
    return 2*x-1

def forward_deriv_f(x,delta):
    return (f(x+delta)-f(x))/delta

x=1
true_value=deriv_f(x)
delta_unicode='\u03b4'
for i in range(-2,-25,-1):
    delta=5**i
    num = forward_deriv_f(x,delta)
    error = abs(num - true_value)
    print(delta_unicode + " = %.2E , Numerical derivative = %.16f , Error = %.2E" %(delta,num,error))
