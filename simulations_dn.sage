from sage.stats.distributions.discrete_gaussian_integer import         DiscreteGaussianDistributionIntegerSampler as DGaussZ

# import Python packages
from mpmath import mpf, sqrt, jtheta, exp
from math import comb, pi
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

def theta_dn(n,s):
    '''
    Compute the theta series of Dn

    :param n: lattice dimension
    :param s: sampling width

    :returns: a real value
    '''
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta_series = mpf('0')
    for m in range(1,n/2+1):
        theta_series += comb(n,2*m)*(theta3**(n-2*m))*(theta2**(2*m))
    return theta_series

def coset_prob_dn(n, s):
    '''
    Compute the probability that a vector lies in a 
    coset of Dn represented by a codeword of weight w=1,...,n/2

    :param n: lattice dimension
    :param s: sampling width

    :returns: a list of real values
    '''
    table=[]
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    for m in range(1,n/2+1):
        p_2m = comb(n,2*m)*(theta3**(n-2*m))*(theta2**(2*m))/theta_dn(n, s)
        table.append(p_2m)
    return table

def sample_dn(n, s, indices, weights_prob):
    ''''
    Sample a lattice vector from Dn 

    :param n: lattice dimension
    :param s: sampling width
    :param indices: a list of vector indices
    :param weights_prob: a list of the probabilities of sampling in each coset, computed offline by coset_prob_dn
    
    :returns: a lattice vector
    '''
    m = random.choices([i for i in range(1,n/2+1)], weights_prob)
    subset = random.sample(indices, k=2*m[0])
    x = np.empty(n)
    for j in range(len(subset)):
        x[subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)()
    not_subset = list(set(indices)-set(subset))
    for j in range(len(not_subset)):
        x[not_subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))()
    return x

def espitau_dn(n, s):
    ''''
    Sample a lattice vector from Dn 

    :param n: lattice dimension
    :param s: sampling width
    
    :returns: a lattice vector
    '''
    x = np.zeros(n)
    x[0] += 1
    while (sum(x) % 2) != 0:
        x = np.array([DGaussZ(sigma=float(s)/(float(sqrt(2*pi))))() for i in range(n)])
    return x

'''
Compute the time it takes to sample 100000 times with both samplers and compare. 
Sampling width s is computed with epsilon = 2**(-36) and we take the example of n = 8
'''

# offline stage
indices = [i for i in range(8)]
weights = coset_prob_dn(8, 2.97)

# online stage
start_time1 = time.perf_counter ()
for i in range(100000):
    sample_dn(8, 2.97, indices, weights)
end_time1 = time.perf_counter ()
print(end_time1 - start_time1, "seconds")

time1 = end_time1 - start_time1

start_time2 = time.perf_counter ()
for i in range(100000):
    espitau_dn(8,2.97)
end_time2 = time.perf_counter ()
print(end_time2 - start_time2, "seconds")

time2 = end_time2 - start_time2

print("Our sampler is", time2/time1, "times as fast as the sampler by Espitau et. al.")