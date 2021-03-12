#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 17:02:45 2021

@author: gumenghan
"""

import math
import numpy as np
import statistics as ss

def monteCarlo(t,n):
    # n is the number of trajactories #
    X = np.random.normal(0,math.sqrt(t),n)
    res = []
    for item in X:
        tmp = item / (max(1,math.pow(item,2)))
        res.append(tmp)
    return res , ss.mean(res)


output = monteCarlo(1.24,100)

# Mean at time 1.24
print(output[1])
