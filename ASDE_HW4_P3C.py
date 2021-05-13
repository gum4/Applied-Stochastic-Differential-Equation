#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 05:58:16 2021

@author: gumenghan
"""

import math
import numpy as np
import statistics as ss
import time
import matplotlib.pyplot as plt



# Simulation for part b
t_end = 1
x_init = 1
alpha = 0.1
beta = 0.3

Num = 10
def dw(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(loc=0.0, scale=np.sqrt(delta))

def Euler_mara_MC(N): # N is the number of time steps
    dt = float(t_end) / N
    ts = np.arange(0, t_end + dt, dt)
    xs = np.zeros(N + 1)

    
    xs[0] = x_init

    
    
    res1= []
    for itr in range(Num):
        for i in range(1, ts.size):
            t =  (i - 1) * dt
            x = xs[i - 1]

            xs[i] = x +alpha * x * dt + beta * dw(dt)
        res1.append(xs)
        plt.plot(ts, xs,'r')
        #plt.plot(ts, ys,'g')
        

    return res1


Euler_mara_MC(1000)