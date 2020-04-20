approx= '\u2248'

def binding_energy(A,Z):
    a1=15.67
    a2=17.23
    a3=0.75
    a4=93.2
    a5=0
    if(A%2==0):
        if(Z%2==0):
            a5=12
        else:
            a5=-12
    return a1*A - a2*(A**(2/3)) - a3*Z*Z*(A**(-1/3))- a4*(A-2*Z)*(A-2*Z)/A + a5*(A**-0.5)


def binding_energy_per_nucleon(A,Z):
    return binding_energy(A,Z)/A

def stable_nucleus(Z):
    maxim=binding_energy_per_nucleon(Z,Z)
    A=Z
    for i in range(Z+1,3*Z+1):
        cur=binding_energy_per_nucleon(i,Z)
        if(cur>maxim):
            A=i
            maxim=cur
    return A,maxim

A = int(input('A = '))
Z = int(input('Z = '))
print("Binding Energy "+ approx +" %.3f MeV " % binding_energy(A,Z))
print("Binding Energy per nucleon "+ approx +" %.3f MeV\n " % binding_energy_per_nucleon(A,Z))
A_opt,maxim=stable_nucleus(Z)
print("Maximum Binding Energy per nucleon (Z=%i) " %Z + approx +" %.3f MeV at A = %i\n " % (maxim,A_opt))

maxim=0
for i in range(1,101):
    A,cur=stable_nucleus(i)
    if(maxim<cur):
        Z_opt=i
        maxim=cur
    print("Z=%i \t A=%i \t  B/A = %.3f MeV" % (i,A,cur))
    
print("Maximum value of binding energy per nucleon occurs at Z = %i" % Z_opt)
