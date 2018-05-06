# -*- coding: utf-8 -*-
"""
Created on Sun May  6 12:18:00 2018

@author: chenyushao
"""

import scipy
from scipy import cos
from scipy import sin
from scipy import sqrt
from scipy import linalg as LA
from scipy import array
from scipy import pi
from scipy import arange
import matplotlib.pyplot as plt
# constants
global hbar, m0, T0, k0
a0 = 5.64
T0 = array([1, 1 ,1])
T0 = T0 * pi / 4
k0 = array([1, 0, 0])
k0 = k0
hbar = 1
m0 = 1

# 15 vectors
# <0, 0, 0>
v0 = array([0, 0, 0])
# <1, 1, 1>
v1 = array([1, 1, 1])
v2 = array([-1, 1, 1])
v3 = array([1, -1, 1])
v4 = array([1, 1, -1])
v5 = array([-1, -1, 1])
v6 = array([1, -1, -1])
v7 = array([-1, 1, -1])
v8 = array([-1, -1, -1])
# <2, 0, 0>
v9 = array([2, 0, 0])
v10 = array([0, 2, 0])
v11 = array([0, 0, 2])
v12 = array([-2, 0, 0])
v13 = array([0, -2, 0])
v14 = array([0, 0, -2])

vectors = [v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14]

def eq(a, b):
    if abs(a - b) < 1e-3:
        return True
    else:
        return False
    
def H(G1, G2, k, T):
    
    q = (G1 - G2).dot(G1 - G2)
    
    if scipy.array_equal(G1, G2):
        h0 = (hbar * hbar / (2 * m0)) * 2 * 4.72 * (G2 + k).dot(G2 + k)
    else:
        h0 = 0
    hS = Vfs(q) * cos((G2 - G1).dot(T))
    hA = sqrt(-1) * (Vfa(q) * sin((G2 - G1).dot(T)))
    return h0 + hS + hA

def Vfs(q):
    value = 0
    if eq(q, (3)):
        value = -3.13
    elif eq(q, (8)):
        value = 0.14
    elif eq(q, (11)):
        value = 0.82
    return value

def Vfa(q):
    value = 0
    if eq(q, (3)):
        value = 0.95
    elif eq(q, (4)):
        value = 0.68
    elif eq(q, (11)):
        value = 0.14
    return value

def energys(k):
    k = k * k0
    Hamiltonian = [[0 + 0*sqrt(-1) for i1 in range(15)] for i2 in range(15)]
    for i1 in range(15):
        for i2 in range(15):
            Hamiltonian[i1][i2] = H(vectors[i1], vectors[i2], k, T0) 
    Hamiltonian = array(Hamiltonian)
    eigs = LA.eigvals(Hamiltonian)
    eigs = scipy.real(eigs)
    Eks = [eigs[i] for i in range(15)]
    Eks.sort()
    return Eks

print(eq(0, 0.0001))

ks = [k for k in arange(0, 1, 0.01)]
krange = range(len(ks))
Ek=[[] for i in range(15)]
for i in range(15):
    Ek[i] = [(energys(k))[i] for k in ks]
    plt.plot(krange, Ek[i], '-')
plt.show()