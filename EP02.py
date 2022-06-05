from cmath import pi
from string import digits
from tkinter import Y
import numpy as np
import math

# W6 W8 W10
W = [[0],
     [2],
     [1],
     [0.8888888888888888888888889, 0.5555555555555555555555556],
     [0.6521451548625461426269360, 0.3478548451374538573730639],
     [],
     [0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961],
     [],
     [0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314],
     [],
     [0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688]]
# X6 X8 X10
X = [[0],
     [0],
     [0.5773502691896257645091487],
     [0, 0.7745966692414833770358531],
     [0.3399810435848562648026658, 0.8611363115940525752239465],
     [],
     [0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016],
     [],
     [0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609],
     [],
     [0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640]]

# Integral Numérica de Gauss: f(x)dx em [-1,1]
def integral_gauss_1A(f, n):
    somatorio = 0
    # Quantidade par de nós:
    if (n%2 == 0):
        for j in range(0, round(n/2)):
            somatorio += (W[n][j]) * (f(-X[n][j]) + f(X[n][j]))
    # Quantidade ímpar de nós:
    else:
        for j in range(0, round(n/2)):
            if (j == 0):
                somatorio += W[n][j] * f(X[n][j])
            else:
                somatorio += (W[n][j]) * (f(-X[n][j]) + f(X[n][j]))
    return somatorio

# Integral Numérica de Gauss: f(x)dx em [a,b]
def integral_gauss_1B(f, a, b, n):
    def g(u):
        return (f((b+a)/2 + (b-a)*u/2) * (b-a)/2)
    return integral_gauss_1A(g, n)

# Integral Numérica de Gauss: f(x,y)dydx em [-1,1] e [-1,1]
def integral_gauss_2A(f, n): 
    somatorio = 0
    # n par:
    if (n%2 == 0):
        for i in range(0, round(n/2)):
            for j in range(0, round(n/2)):
                somatorio += (W[n][i]*W[n][j] * (f(-X[n][i],-X[n][j]) + f(X[n][i],X[n][j])))
    # n ímpar:
    else:
        for i in range(0, round(n/2)):
            for j in range(0, round(n/2)):
                if (j == 0):
                    somatorio += W[n][i]*W[n][j] * f(X[n][i], X[n][j])
                else:
                    somatorio += (W[n][i]*W[n][j] * (f(-X[n][i],-X[n][j]) + f(X[n][i],X[n][j])))
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
    somatorio = 0
    # n par:
    if (n%2 == 0):
        for i in range(0, round(n/2)):
            for j in range(0, round(n/2)):
                for k in range(0, round(n/2)):
                    somatorio += (W[n][i]*W[n][j]*W[n][k] * (f(-X[n][i],-X[n][j],-X[n][k]) + f(X[n][i],X[n][j],X[n][k])))
    # n ímpar:
    else:
        for i in range(0, round(n/2)):
            for j in range(0, round(n/2)):
                for k in range(0, round(n/2)):
                    if (j == 0):
                        somatorio += (W[n][i]*W[n][j]*W[n][k] * f(X[n][i],X[n][j],X[n][k]))
                    else:
                        somatorio += (W[n][i]*W[n][j]*W[n][k] * (f(-X[n][i],-X[n][j],-X[n][k]) + f(X[n][i],X[n][j],X[n][k])))
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
    def f(x,y,z):
        return 1
    resultados = np.array([])
    n = [6, 8, 10]
    for i in n:
        resultados = np.append(integral_gauss_3B(f,0,1,0,1,0,1,i))
    print("\nExemplo 1A - Área do cubo de aresta 1:\n")
    for i in range(3):
        print('\tn = {}\tResultado = {:0.22f}'.format(n[i], resultados[i]))

def main():
    exemplo_1A()
main()
