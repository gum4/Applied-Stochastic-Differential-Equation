#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:09:31 2021

@author: gumenghan
"""
import numpy as np
import matplotlib as plt
import statistics as ss
# Monte Carlo Simulation


# num is the number of iterations (time steps) #

def monteCarlo(a,b,sigmaZ,sigma0,mu0,num):
    # n is the number of trajactories #
    n = 10 ** 3
    X0 = np.random.normal(mu0,sigma0,1)
    x0 , res = X0[0] , []
    for i in range(n): 
        Random = np.random.normal(0,sigmaZ,num)
        xn=x0
        for item in Random:
            C = a * (xn - b)
            xn += C + item
        res.append(xn)
    return res


# When a=-1 as time step increases the mean goes to constant b #
a, b, sigmaZ, sigma0, mu0 = -1.2, 3, 0.2, 0.1, 1
E, Var = [], []
for num in range (1,10 ** 2):
    tmp = monteCarlo(a, b, sigmaZ, sigma0, mu0, num)
    E.append(ss.mean(tmp))
    Var.append(ss.variance(tmp))

print (E)

# Calculation of the covariable between 2 time steps #$
i, j= 36, 77
n1 = monteCarlo(a, b, sigmaZ, sigma0, mu0, i)
n2 = monteCarlo(a, b, sigmaZ, sigma0, mu0, j)
sum = 0

for i in range(0, len(n1)):
    sum += ((n1[i] - ss.mean(n1)) * (n2[i] - ss.mean(n2)))

print(sum/(len(n1)-1))

            

            

