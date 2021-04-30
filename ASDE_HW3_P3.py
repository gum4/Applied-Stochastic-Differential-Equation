#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 05:23:33 2021

@author: gumenghan
"""


import numpy as np
import statistics as ss
import matplotlib.pyplot as plt

# Display N runs
Num = 20

t_init = 0 #t_init must be 0
t_end  = 1.5

theta = 0.3 # mean in GBM
Sigma = 0.05 # SDV in GBM
y_init = 10 # initial value at t_init

def Tao(t):
    return 0.5*t

def f(x):
    return 0.02 * x

def theta(t):
    return 2 * t

def sigma(t):
    return t

# W(t_n+1) - W(t_n)
def dW(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(loc=0.0, scale=np.sqrt(delta))

def Euler_mara(N):
    t1 = t_init + Tao(t_end - t_init)
    dt = float(t_end - t_init) / N
    ts = np.arange(t_init, t_end + dt, dt)
    ys = np.zeros(N + 1)

    ys[0] = y_init
    res = []
    for itr in range(Num):
        for i in range(1, ts.size):
            t = t_init + (i - 1) * dt
            y = ys[i - 1]
            ys[i] = y + theta (t) * y * dt + sigma (t) * y * dW(dt)
        cur = ys[ts.size-1]
        pre = ys[int(ts.size * t1 / (t_end - t_init)) - 1]
        share = f (pre)
        res.append(share * (cur - y_init) )
        
    return res


res = Euler_mara(1000)
print(ss.mean(res))
print(ss.variance(res))