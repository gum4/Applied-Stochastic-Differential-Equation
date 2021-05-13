#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 01:06:49 2021

@author: gumenghan
"""

import math
import numpy as np
import statistics as ss
import time
import matplotlib.pyplot as plt

Num=10
t_end = 1
epsilon = 0.01
sigma = 0.3
theta = 2
x_init=1
y_init=1

# W(t_n+1) - W(t_n)
def dW(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(0, np.sqrt(delta),2)

def dw(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(loc=0.0, scale=np.sqrt(delta))

def PHI_dx(x,y): #partial derivative of PHI(x,y) on x
    return 2*math.pow(x,1)*math.pow(y,2)

def PHI_dy(x,y): #partial derivative of PHI(x,y) on y
    return 2*math.pow(x,2)*math.pow(y,1)

def PHI_prime_x(x): #based on theta=2
    return math.log(x) - math.log(math.pi)/2

def PHI_prime_dx(x):#based on theta=2
    return 1.0/x

# Store all X(t) instead of X(t) at end_time in two systems
def Euler_mara_system1(N): # N is the number of time steps
    dt = float(t_end) / N
    ts = np.arange(0, t_end + dt, dt)
    xs = np.zeros(N + 1)
    ys = np.zeros(N + 1)
    
    xs[0] = x_init
    ys[0] = y_init
    
    
    res1= []
    for itr in range(Num):
        for i in range(1, ts.size):
            t =  (i - 1) * dt
            x = xs[i - 1]
            y = ys[i - 1]
            DW = dW(dt)
            xs[i] = x - PHI_dx(x,y) * dt + math.sqrt(2 * sigma) * DW[0]
            ys[i] = y - PHI_dy(x,y) * dt / epsilon + math.sqrt(2 * theta / epsilon) * DW[1]
        res1.append(xs)
        plt.plot(ts, xs,'r')
        #plt.plot(ts, ys,'g')
        

    return res1



def Euler_mara_system2(N): # N is the number of time steps
    dt = float(t_end) / N
    ts = np.arange(0, t_end + dt, dt)
    xs = np.zeros(N + 1)
    
    xs[0] = x_init
    
    res1 = []
    for itr in range(Num):
        for i in range(1, ts.size):
            t =  (i - 1) * dt
            x = xs[i - 1]
            xs[i] = x - PHI_prime_dx(x) * dt + math.sqrt(2 * sigma) *dw(dt)
        res1.append(xs)
        plt.plot(ts, xs,'g')
        #plt.plot(ts, ys,'g')
        
    plt.show()
    return res1


plt.xlabel("time (s)")
h = plt.ylabel("X(t)")
h.set_rotation(0)
t11=time.time()
E1 = Euler_mara_system1(1000)
print(time.time()-t11)
t12=time.time()
E2 = Euler_mara_system2(1000)
print(time.time()-t12)
E1, E2 = np.array(E1), np.array(E2)

plt.show()


# Compare the mean and covariance of system <7> and <8>
print(E1.mean(axis=0))
print(E2.mean(axis=0))
print(np.cov(E1,rowvar=False))
print(np.cov(E2,rowvar=False))

