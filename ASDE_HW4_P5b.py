import math
import numpy as np
import statistics as ss
import time
t_end=1
X_init=0
V_init=0
Num=100
mu=1

def PHI_prime(x):
    return  math.exp(x)   -x +1

def dW(delta):
    #Sample a random number
    # standard deviation is sqrt(delta), delta is time interval
    return np.random.normal(loc=0.0, scale=np.sqrt(delta))

# N is the number of grid points, fix starting time and end time, increase N to reduce delta_t
def Euler_mara(N,theta,a):
    dt = float(t_end) / N
    ts = np.arange(0, t_end + dt, dt)
    X = np.zeros(N + 1)
    V = np.zeros(N + 1)
    tau = np.zeros(Num)
    X[0] = X_init
    V[0] = V_init
    
    for itr in range(Num):
        flag=False
        for i in range(1, ts.size):
            t = (i - 1) * dt
            x = X[i - 1]
            v = V[i - 1]
            X[i] = x + math.pow(theta,-1) * v * dt
            V[i] = v - math.pow(theta,-2) * v * dt - (1.0/theta) * mu * PHI_prime(x) * dt + (math.sqrt(2) / theta ) * dW(dt)
            
        # COmpute first passage time
        for i in range (1,ts.size):
            if X[i]<=a:
                tau[itr] = 0.5 * dt * (i + i-1)
                flag = True
                break
        if not flag:
            tau[itr]=t_end

        
    return tau

t1 = time.time()
res=Euler_mara(10000, 1, -10)
print (time.time()-t1) # CPU time for running Euler Marayama method
print(ss.mean(res), ss.stdev(res))
