#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 18:46:35 2021

@author: gumenghan
"""

import math
import numpy as np
import statistics as ss
import time
import matplotlib.pyplot as plt

t_end = 1
Num = 100  #simulate Num paths

def dW(delta,D):
    #Sample D random numbers
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(0, np.sqrt(delta),D)


def Hypersperical_monteCarlo(N,D,rho,start): # start is r0 D is dimension, N is time steps
    X_init=np.full(D,start)
    
    dt = float(t_end) / N
    ts = np.arange(0, t_end + dt, dt)
    
    tau = np.zeros(Num) #Expected time to reach outer boundary
    
    
   
    for itr in range(Num):
        X = np.zeros((D ,  N + 1))
        X[:,0] = X_init
        flag=False # checks whether reach a boundary before t_end
        for i in range(1, ts.size):
            dw = dW(dt, D)
            t = (i - 1) * dt
            x = X[:,i - 1]
            for j in range(D):
                # increment in D-dimensional brownian motion
                X[j,i] = x[j] + dw[j]
            
        # Compute Expected time to reach outer boundary
        for i in range (1,ts.size):
            Sum=0
            for j in range(D):
                Sum += math.pow(X[j,i] - X_init[j], 2)

            if Sum >= math.pow((rho -start),2):
                tau[itr] = 0.5 * dt * (i + i-1)
                flag = True
                break
        if not flag:
            # if not found, first passage time is t_end
            tau[itr]=t_end
        
        
    return tau

res1_1 = Hypersperical_monteCarlo(1000, 1, 2, 1.5) #ro=1.5 outer boundary is 2
res1_2 = Hypersperical_monteCarlo(1000, 1, 1, 1.5) #ro=1.5 inner boundary is 1

print (ss.mean(res1_1), ss.mean(res1_2))

res2_1 = Hypersperical_monteCarlo(1000, 2, 2, 1.5) #ro=1.5 outer boundary is 2
res2_2 = Hypersperical_monteCarlo(1000, 2, 1, 1.5) #ro=1.5 inner boundary is 1

print (ss.mean(res2_1), ss.mean(res2_2))

res10_1 = Hypersperical_monteCarlo(1000, 10, 2, 1.5) #ro=1.5 outer boundary is 2
res10_2 = Hypersperical_monteCarlo(1000, 10, 1, 1.5) #ro=1.5 inner boundary is 1

print (ss.mean(res10_1), ss.mean(res10_2))