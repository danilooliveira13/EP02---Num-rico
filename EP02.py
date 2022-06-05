from cmath import pi
from string import digits
from tkinter import Y
import numpy as np
import math

from sympy import root

# W6 W8 W10
W_d = [[],
     [],
     [],
     [],
     [],
     [],
     [0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961],
     [],
     [0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314],
     [],
     [0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688]]
# X6 X8 X10
X_d = [[],
     [],
     [],
     [],
     [],
     [],
     [0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016],
     [],
     [0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609],
     [],
     [0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640]]

# Integral Numérica de Gauss: f(x)dx em [-1,1]
def integral_gauss_1A(f, n):
    X = np.array(X_d[n])
    W = np.array(W_d[n])
    X = np.concatenate((-X[::-1], X))
    W = np.concatenate((W[::-1], W))
    somatorio = 0
    for i in range(n):
        somatorio += W[i] * f(X[i])    
    return somatorio

# Integral Numérica de Gauss: f(x)dx em [a,b]
def integral_gauss_1B(f, a, b, n):
    def g(u):
        return (f((b+a)/2 + (b-a)*u/2) * (b-a)/2)
    return integral_gauss_1A(g, n)

# Integral Numérica de Gauss: f(x,y)dydx em [-1,1] e [-1,1]
def integral_gauss_2A(f, n):
    X = np.array(X_d[n])
    W = np.array(W_d[n])
    X = np.concatenate((-X[::-1], X))
    W = np.concatenate((W[::-1], W))
    somatorio = 0
    for i in range(n):
        for j in range(n):
            somatorio += (W[i]*W[j] * f(X[i],X[j]))
    return somatorio

# Integral Numérica de Gauss: f(x,y)dydx em [a,b] e [c(x),d(x)]
def integral_gauss_2B(f, a, b, c, d , n):
    def g1(x,v):
        return (f(x,(d(x)+c(x))/2 +(d(x)-c(x))*v/2) * (d(x)-c(x))/2)
    def g2(u,v):
        return (g1((b+a)/2 + (b-a)*u/2,v) * (b-a)/2)   
    return integral_gauss_2A(g2, n)

# Integral Numérica de Gauss: f(x,y,z)dzdydx em [-1,1], [-1,1] e [-1,1]
def integral_gauss_3A(f, n):
    X = np.array(X_d[n])
    W = np.array(W_d[n])
    X = np.concatenate((-X[::-1], X))
    W = np.concatenate((W[::-1], W))
    somatorio = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                somatorio += (W[i]*W[j]*W[k] * f(X[i],X[j],X[k]))
    return somatorio 

# Integral Numérica de Gauss: f(x,y,z)dzdydx em [a,b], [c(x),d(x)] e [e(x,y),f(x,y)]
def integral_gauss_3B(f1, a, b, c, d, e, f, n):
    def g1(x,y,w):
        return (f1(x,y,(f(x,y)+e(x,y))/2 + (f(x,y)-e(x,y))*w/2) * (f(x,y)-e(x,y))/2)
    def g2(x,v,w):
        return (g1(x,(d(x)+c(x))/2 +(d(x)-c(x))*v/2,w) * (d(x)-c(x))/2)
    def g3(u,v,w):
        return (g2((b+a)/2 + (b-a)*u/2,v,w) * (b-a)/2)
    return integral_gauss_3A(g3, n)

# Exemplo 1A: área do cubo
def exemplo_1A():
    def f1(x,y,z):
        return 1
    def c(x):
        return 0
    def d(x):
        return 1
    def e(x,y):
        return 0
    def f(x,y):
        return 1
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_3B(f1,0,1,c,d,e,f,i))
    print("\nExemplo 1A - Área do cubo de aresta 1:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_1B():
    def f1(x,y,z):
        return 1
    def c(x):
        return 0
    def d(x):
        return 1-x
    def e(x,y):
        return 0
    def f(x,y):
        return 1-x-y
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_3B(f1,0,1,c,d,e,f,i))
    print("\nExemplo 1B - Área do tetraedro:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_2A():
    def f1(x,y):
        return 1
    def c(x):
        return 0
    def d(x):
        return 1-x**2
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_2B(f1,0,1,c,d,i))
    print("\nExemplo 2A - Área no primeiro quadrante, limitada pelos eixos e por y = 1 - x^2:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_2B():
    def f1(y,x):
        return 1
    def c(y):
        return 0
    def d(y):
        return root(1-y,2)
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_2B(f1,0,1,c,d,i))
    print("\nExemplo 2B - Área no primeiro quadrante, limitada pelos eixos e por y = 1 - x^2:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_3A():
    def f1(x,y):
        fx = -(np.e**(y/x))*(y/x**2)
        fy = (np.e**(y/x))/x
        func = root(fx**2  + fy**2 +1, 2)
        return func
    def c(x):
        return x**3
    def d(x):
        return x**2
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_2B(f1,0.1,0.5,c,d,i))
    print("\nExemplo 3A - Área da superfíe descrita por z = e^(y/x):")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_3B():
    def f1(x,y,z):
        return 1
    def c(x):
        return x**3
    def d(x):
        return x**2
    def e(x,y):
        return 0
    def f(x,y):
        return np.e**(y/x)
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_3B(f1,0.1,0.5,c,d,e,f,i))
    print("\nExemplo 3B - Volume da superfíe descrita por z = e^(y/x):")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_4A():
    def f1(x,y):
        return 2*pi*y
    def c(x):
        return 0
    def d(x):
        return root(1/4 - x**2, 2)
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_2B(f1,0,1/4,c,d,i))
    print("\nExemplo 4A - Volume da calota esférica de altura 1/4 da altura da esfera de raio 1:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def exemplo_4B():
    def f1(y,x):
        return 2*pi*x
    def c(y):
        return 0
    def d(y):
        return np.e**(-y**2)
    resultados = []
    n = [6, 8, 10]
    for i in n:
        resultados.append(integral_gauss_2B(f1,-1,1,c,d,i))
    print("\nExemplo 4A - Volume do sólido de revolução:")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def main():
    print('\nExercício Programa 2: INTEGRAÇÃO NUMÉRICA DE GAUSS')
    print('Autores:\n\tRenato Naves Fleury\n\tDanilo Oliveira da Silva')
    exemplo_1A()
    exemplo_1B()
    exemplo_2A()
    exemplo_2B()
    exemplo_3A()
    exemplo_3B()
    exemplo_4A()
    exemplo_4B()
    print('\n')
main()
