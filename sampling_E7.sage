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
    
def theta_E7(s):
    '''
    Compute the theta series of E7

    :param s: sampling width

    :returns: a real value
    '''
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta_E7 = theta3**7+7*(theta3**4)*(theta2**3)+7*(theta3**3)*(theta2**4)+theta2**7
    return theta_E7

def coset_prob_E7(s):
    '''
    Compute the probability that a vector lies in a 
    coset of E7 represented by a codeword of weight w=0,3,4,7

    :param s: sampling width

    :returns: a list of real values
    '''
    table=[]
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta_series = theta_E7(s)
    coset_0 = theta3**7/theta_series
    coset_3 = 7*(theta3**4)*(theta2**3)/theta_series
    coset_4 = 7*(theta3**3)*(theta2**4)/theta_series
    coset_7 = theta2**7/theta_series
    table.extend([coset_0, coset_3, coset_4, coset_7])
    return table

# Precompute the code associated with E7 and the weight 3 and 4 codewords (i.e. the non trivial codewords)

C = codes.HammingCode(GF(2),3)
weight_3 = [[1, 0, 0, 0, 0, 1, 1], [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0],
            [1, 1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0], [0, 1, 0, 1, 0, 1, 0],
            [0, 0, 1, 1, 0, 0, 1]]

weight_4 = [[1, 1, 0, 0, 1, 1, 0], [1, 0, 1, 0, 1, 0, 1], [0, 1, 1, 0, 0, 1, 1],
            [0, 0, 0, 1, 1, 1, 1], [1, 1, 0, 1, 0, 0, 1], [1, 0, 1, 1, 0, 1, 0],
            [0, 1, 1, 1, 1, 0, 0]]

def sample_E7(s, weights_prob):
    ''''
    Sample a lattice vector from E7 

    :param s: sampling width
    :param weights_prob: a list of the probabilities of sampling in each coset, computed offline by coset_prob_E7
    
    :returns: a lattice vector
    '''
    w = random.choices([0,3,4,7], weights=weights_prob)
    if w[0] == 3:
        c = random.sample(weight_3, k=1)
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))), c=c[0][i]/2)() for i in range(7)])
    elif w[0] == 4:
        c = random.sample(weight_4, k=1)
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))), c=c[0][i]/2)() for i in range(7)])
    elif w[0] == 0: 
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))() for i in range(7)])
    elif w[0] == 7:
        x = np.array([2*DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)() for i in range(7)])
    return x

'''
Compute the time it takes to sample 100000 times with both samplers and compare. 
Sampling width s is computed with epsilon = 2**(-36)
'''

start_time1 = time.perf_counter ()
for i in range(100000):
    sample_E7(1.75, coset_prob_E7(1.75))
end_time1 = time.perf_counter ()
print(end_time1 - start_time1, "seconds")