import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#-------- Functions -------#

   
def NLSE_t(t, f, beta, K, G, w0):
    """ non-linear ODE. f is usually complex.
        n is the index of the wave guides
    """   
   
    N = len(f)
    # random numbers from -1 to 1
    rand = (random(N)-0.5)*2.
    a =  beta * (f * np.conj(f) * f) -G*rand*f
    for i in range(1,N-1):
        a[i] += K * (f[i-1]  - 2.0 * f[i]  + f[i+1]) #
    a[0]     += K * (f[N-1]  - 2.0 * f[0] + f[1])   
    a[N-1]   += K * (f[0]  - 2.0 * f[N-1]  + f[N-2])  
    return a / (-1.* w0 * 1j)



# ----  Create arrays ------#

# ------ Parameters --------#

K = 1. # hopping amplitude between waveguides
beta = 1. # Non-linearity (factor 3 in ODE)
G = 0.#1e-5,2.5e-5,5e-5,1e-4,2.5e-4,5e-4,1e-3
N = 384 # number of waveguide
w0 = 1. # scaling 
T = 1592 # Length of system
phi = 0.5 #1. or 0.5 in paper
alpha = 0.001   #0.05 in paper
P =  32#30
time_t=np.arange(0,l)
p =4 #100
t0 = 0
methods = ['DBF','RK45'] #numerical method for integration
rtols = [1e-9, 1e-8] 
use_both = False

# compute other properties
l = T * 2. * np.pi  # length of system
A = 2. * phi
#q and Q should be between 0 and pi
assert p >= 0
assert p <= N//2
assert P >= 0
assert P <= N//2
Q = (2. * np.pi * P) / N
q = (2. * np.pi * p) / N
w = np.sqrt(w0**2 + 4 * K * np.sin(q/2.)**2 + G )
dtmax = 2 * np.pi / max(w, w0) / 50.


# set initial values
nn = np.arange(N)
a0 = np.zeros(N, dtype=complex)
a0[:] = (A + alpha * np.exp(1j*Q*nn)) * np.exp(1j*(q *nn - w * t0))





if 1:
    # generate data

    sol = solve_ivp(NLSE_t, [0,l], a0, method=methods[1], t_eval=time_t, rtol=rtols[1], args=(beta, K, G, w0))
    np.savez("BDF_U1.npz", t=sol.t, y=sol.y)
    t = sol.t
    y = sol.y


