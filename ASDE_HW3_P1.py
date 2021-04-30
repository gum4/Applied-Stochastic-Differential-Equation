#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 13:52:17 2021

@author: gumenghan
"""



import numpy as np
import statistics as ss
import matplotlib.pyplot as plt

# Display N runs
Num = 20

t_init = 0
t_end  = 10

theta = 0.3 # mean in GBM
Sigma = 0.05 # SDV in GBM
y_init = 10 # initial value at t_init



# W(t_n+1) - W(t_n)
def dW(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(loc=0.0, scale=np.sqrt(delta))

# N is the number of grid points, fix starting time and end time, increase N to reduce delta_t
def Euler_mara(N):
    dt = float(t_end - t_init) / N
    ts = np.arange(t_init, t_end + dt, dt)
    ys = np.zeros(N + 1)

    ys[0] = y_init
    plt.xlabel("time (s)")
    h = plt.ylabel("y")
    h.set_rotation(0)
    res, f_res = [], []
    for itr in range(Num):
        for i in range(1, ts.size):
            t = t_init + (i - 1) * dt
            y = ys[i - 1]
            ys[i] = y + theta * y * dt + Sigma * y * dW(dt)
        cur = ys[ts.size-1]
        res.append(cur)
        f_res.append(2*cur)
        plt.plot(ts, ys)
        
    plt.show()
    return res, f_res

def Trapezoidal(N):
    dt = float(t_end - t_init) / N
    ts = np.arange(t_init, t_end + dt, dt)
    ys = np.zeros(N + 1)

    ys[0] = y_init
    plt.xlabel("time (s)")
    h = plt.ylabel("y")
    h.set_rotation(0)
    res = []
    for itr in range(Num):
        for i in range(1, ts.size):
            t = t_init + (i - 1) * dt
            y = ys[i - 1]
            y_hat = y + theta * y * dt + Sigma * y * dW(dt)
            ys[i] = y + 0.5 * (theta * y + theta * y_hat) * dt + 0.5 * (Sigma * y +Sigma * y_hat) * dW(dt)
        cur = ys[ts.size-1]
        res.append(cur)
        plt.plot(ts, ys)
        
    plt.show()
    return res


target = y_init * np.exp (theta * (t_end - t_init))

# Euler Marayama is strongly and weakly convergent
Monte1 = Euler_mara(10)
diff1_strong = ss.mean(Monte1[0]) - target
diff1_weak = ss.mean(Monte1[1]) - 2*target
print(diff1_strong, diff1_weak)
Monte2 = Euler_mara(1000)
diff2_strong = ss.mean(Monte2[0]) - target
diff2_weak = ss.mean(Monte2[1]) - 2*target
print(diff2_strong, diff2_weak)
Monte3 = Euler_mara(20000)
diff3_strong = ss.mean(Monte3[0]) - target
diff3_weak = ss.mean(Monte3[1]) - 2*target
print(diff3_strong, diff3_weak)

# Trapezoidal method is not convergent
Trape1 = Trapezoidal(10)
Diff1 = ss.mean(Trape1) - target
Trape2 = Trapezoidal(1000)
Diff2 = ss.mean(Trape1) - target
Trape3 = Trapezoidal(20000)
Diff3 = ss.mean(Trape1) - target
print (Diff1, Diff2, Diff3)
