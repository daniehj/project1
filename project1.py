from pylab import *
import time

def tridiagonal(n):
    
    time.clock()

    x = linspace(0,1,n+2)
    
    a = zeros(n+2)
    b = zeros(n+2)
    c = zeros(n+2)
    
    a.fill(-1)
    b.fill(2)
    c.fill(-1)
    
    h = x[2] - x[1]
    
    f = h**2*100*exp(-10*x)
    u = 1. - (1. - exp(-10))*x - exp(-10*x)
    

    f_mark = zeros(n+2)
    v = zeros(n+2)
    
    f_mark[1] = f[1]
    for i in range(2,n):
        b[i] = b[i] - a[i]*c[i-1]/b[i-1]
        f_mark[i] = f[i] - a[i]*f_mark[i-1]/b[i-1]
        
    v[n] = f_mark[i]
    
    while i >= 1:
        v[i] = (f_mark[i] - c[i]*v[i+1])/b[i]
        i-=1
       
    t =  time.clock()
    print 'Time used for n=%.0f , %ss:' % (n,t)


    

##########################################

def tridiagonal_opt(n):
    
    time.clock()
    
    x = linspace(0,1,n+2)
    
    a = zeros(n+2)
    b = zeros(n+2)
    c = zeros(n+2)
    
    a.fill(-1)
    b.fill(2)
    c.fill(-1)
    
    h = x[2] - x[1]
    
    f = h**2*100*exp(-10*x)
    u = 1. - (1. - exp(-10))*x - exp(-10*x)
    
    figure()
    plot(u)
    
    f_mark = zeros(n+2)
    v = zeros(n+2)
    
    f_mark[1] = f[1]
    for i in range(2,n):
        b[i] = b[i] - 1./b[i-1]
        f_mark[i] = f[i] + f_mark[i-1]/b[i-1]
        
    v[n] = f_mark[i]
    
    while i >= 1:
        v[i] = (f_mark[i] + v[i+1])/b[i]
        i-=1
    
    stn = str(n)    
    plot(v,'--r')
    title(r'$n=$' + stn)
    xlabel('')
    ylabel('')
    legend((r'$Analytical$',r'$Numerical$'))
    #savefig('p1n' + stn + '.png')
    close()
    
    err = zeros(n+2)
    
    for i in range(2,n+1):
        err[i] = log10(abs((v[i]-u[i])/u[i]))
        
        
    err = err[2:-1]
    
    return max(err)

def LUdecomp(n):
    
    time.clock()
    
    A = zeros(shape=(n+2,n+2))
    
    x = linspace(0,1,n+2)
    
    a = zeros(n+2)
    b = zeros(n+2)
    c = zeros(n+2)
    
    a.fill(-1)
    b.fill(2)
    c.fill(-1)
    
    h = x[2] - x[1]
    
    f = h**2*100*exp(-10*x)
    u = 1. - (1. - exp(-10))*x - exp(-10*x)
    
    
    for i in range(n+2):
        if (i > 0):
            A[i, i-1] = a[i];
            A[i, i] = b[i];

        if (i < n+1):
            A[i, i+1] = c[i]
    A[0,0] = 2
    v = scipy.linalg.solve(A,f)
    P, L, U = scipy.linalg.lu(A)
    
    t =  time.clock()
    print 'Time used for n=%.0f , %ss:' % (n,t)

tot_err = []
'''
m = [10,100,1000,10000,100000,1000000,10000000]


for n in m:
    tridiagonal(n)




for n in m:
   tot_err.append(tridiagonal_opt(n))

h = 1./(array(m)+1)
figure()
plot(tot_err)
title('Error plot')
xlabel(r'$h**-m$')
ylabel('relative error')
savefig('p1errorplt.png')
show()
'''

'''
Test gives plot and:
C:\FYS\FYS3150>project1.py
Time used for n=10 , 0.000103111233317s:
Time used for n=100 , 0.00113066800672s:
Time used for n=1000 , 0.00535822857272s:
Time used for n=10000 , 0.0354165358023s:
Time used for n=100000 , 0.328818167488s:
Time used for n=1000000 , 3.22548214378s:
Time used for n=10000000 , 32.0908654904s:
'''
m = [10,100,1000,10000]

import scipy.linalg
 
for n in m:
    LUdecomp(n)
 
 '''
 Test gives:
Time used for n=10 , 0.119274808029s:
Time used for n=100 , 0.121753823313s:
Time used for n=1000 , 0.24466122824s:
 '''
