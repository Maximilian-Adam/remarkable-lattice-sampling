from sage.stats.distributions.discrete_gaussian_integer import         DiscreteGaussianDistributionIntegerSampler as DGaussZ

# import Python packages
from mpmath import sqrt, jtheta, exp
from math import pi
import numpy as np
import time
import random 

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

def coset_prob_A2(s):
    '''
    Compute the probability that a vector lies in a coset of A2

    :param s: sampling width

    :returns: a list of real values
    '''
    q = nome(s)
    theta_A2 = jtheta(3,0,q)*jtheta(3,0,q**3)+jtheta(2,0,q)*jtheta(2,0,q**3)
    prob = jtheta(3,0,q)*jtheta(3,0,q**3)/theta_A2
    return [prob, 1-prob]

def sample_A2(s, weights_prob):
    ''''
    Sample a lattice vector from A2 

    :param s: sampling width
    :param weights_prob: a list of the probabilities of sampling in each coset, computed offline by coset_prob_A2
    
    :returns: a lattice vector
    '''
    cosets =[[0,0], [1/2, sqrt(3)/2]]
    t = random.choices(cosets, weights=weights_prob)
    if t == [0,0]:
        x = np.array(DGaussZ(sigma=float(s)/float(sqrt(2*pi)))())
        x = np.append(x, sqrt(3)*DGaussZ(sigma=float(s)/float(sqrt(3*2*pi)))())
    else:
        x = np.array(DGaussZ(float(s)/float(sqrt(2*pi)),c=1/2)())
        np.append(x, sqrt(3)*DGaussZ(sigma=float(s)/float(sqrt(3*2*pi)), c=1/2)())
    return x

'''
Compute the time it takes to sample 100000 times with both samplers and compare. 
Sampling width s is computed with epsilon = 2**(-36)
'''

coset_prob = coset_prob_A2(1.46)

start_time = time.perf_counter ()
for i in range(400000):
    sample_A2(1.46, coset_prob)
end_time = time.perf_counter ()
print(end_time - start_time, "seconds")