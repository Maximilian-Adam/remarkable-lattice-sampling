from sage.stats.distributions.discrete_gaussian_integer import         DiscreteGaussianDistributionIntegerSampler as DGaussZ

# import Python packages
from mpmath import sqrt, jtheta, exp
from math import pi
import numpy as np
import random 
import time

'''
Note:
The Gaussian function with sampling width s is given by rho(x) = exp(-1*pi*||x||**2/s**2). The discrete Gaussian 
ZZ sampler from Sage samples with respect to the Gaussian function exp(-x**2/(2*s**2)), so we rescale the sampling width 
in the input to DGaussZ by 1/sqrt(2*pi). The following code is used to compare samplers for Dn with Espitau et. al. 
'''

def nome(s):
    '''
    Compute the nome q = exp(-1*pi/s**2)

    :param s: sampling width

    :returns: a real value
    '''
    return exp(-1*pi/s**2)

def theta_E8(s):
    '''
    Compute the theta series of E8

    :param s: sampling width

    :returns: a real value
    '''
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta_E8 = theta3**8+14*(theta3**4)*(theta2**4)+theta2**8
    return theta_E8

def coset_prob_E8(s):
    '''
    Compute the probability that a vector lies in a 
    coset of E8 represented by a codeword of weight w=0,4,8

    :param s: sampling width

    :returns: a list of real values
    '''
    table=[]
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta_series = theta_E8(s)
    coset_0 = theta3**8/theta_series
    coset_4 = 14*(theta3**4)*(theta2**4)/theta_series
    coset_8 = theta2**8/theta_series
    table.extend([coset_0, coset_4, coset_8])
    return table

# Precompute the code associated with E8 and the weight 4 codewords (i.e. the non trivial codewords)

C = codes.BinaryReedMullerCode(1,3)

weight_4 = [[0, 1, 0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 0, 1, 1], [1, 1, 0, 0, 1, 1, 0, 0],
            [0, 1, 1, 0, 0, 1, 1, 0], [1, 0, 0, 1, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1, 1, 1], [1, 1, 1, 1, 0, 0, 0, 0],
            [0, 1, 0, 1, 1, 0, 1, 0], [1, 0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 1, 1, 1, 0, 0], [1, 1, 0, 0, 0, 0, 1, 1], 
            [0, 1, 1, 0, 1, 0, 0, 1], [1, 0, 0, 1, 0, 1, 1, 0]]

def sample_E8(s, weights_prob):
    ''''
    Sample a lattice vector from E8 

    :param s: sampling width
    :param weights_prob: a list of the probabilities of sampling in each coset, computed offline by coset_prob_E8
    
    :returns: a lattice vector
    '''
    w = random.choices([0,4,8], weights=weights_prob)
    if w[0] == 4:
        c = random.sample(weight_4, k=1)
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))), c=c[0][i]/2)() for i in range(8)])
    elif w[0] == 0: 
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))() for i in range(8)])
    elif w[0] == 8:
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)() for i in range(8)])
    return x

def rho_vec(v, s):
    ''''
    Compute probability of a vector 

    :param v: vector
    :param s: sampling width
    
    :returns: a real value between 0 and 1
    '''
    v_t = np.transpose(v)
    Sigma = (1/s**2)*np.identity(len(v))
    prod = np.dot(Sigma,v)
    return exp(-1*pi*np.dot(v_t,prod))

def espitau_dn_shift(n, s, t):
    ''''
    Sample a lattice vector from Dn with shifted center

    :param n: lattice dimension
    :param s: sampling width
    :param t: center
    
    :returns: a lattice vector
    '''
    x = np.zeros(n)
    x[0] += 1
    while (sum(x) % 2) != 0:
        x = np.array([DGaussZ(sigma=float(s)/(float(sqrt(2*pi))), c=float(t[i]))() for i in range(n)])
    return x

def espitau_E8(s, s_prime):
    ''''
    Sample a lattice vector from E8 

    :param s: smoothing parameter of D8
    :param s_prime: smoothing parameter of E8
    
    :returns: a lattice vector
    '''
    b = random.sample([0,1], k=1)[0]
    bh = np.array([b*1/2 for i in range(8)])
    x = espitau_dn_shift(8, s,-bh)
    samp = bh+x
    prob = rho_vec(samp, s_prime)/rho_vec(samp, s)
    resp = random.choices(['accept', 'reject'], weights=[prob, 1-prob])[0]
    if resp == 'accept':
        return samp
    else:
        samp = espitau_E8(s, s_prime)
        return samp

'''
Compute the time it takes to sample 100000 times with both samplers and compare. 
Sampling width s is computed with epsilon = 2**(-36)
'''

start_time1 = time.perf_counter ()
for i in range(100000):
    sample_E8(2.2, coset_prob_E8(2.2))
end_time1 = time.perf_counter ()
print(end_time1 - start_time1, "seconds")

time1 = end_time1 - start_time1

start_time2 = time.perf_counter ()
for i in range(100000):
    espitau_E8(2.97, 2.2)
end_time2 = time.perf_counter ()
print(end_time2 - start_time2, "seconds")

time2 = end_time2 - start_time2

print("Our sampler is", time2/time1, "times as fast as the sampler by Espitau et. al.")