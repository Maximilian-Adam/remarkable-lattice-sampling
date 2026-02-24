from sage.stats.distributions.discrete_gaussian_integer import         DiscreteGaussianDistributionIntegerSampler as DGaussZ

# import Python packages
from mpmath import mpf, sqrt, jtheta, exp
from math import comb, pi, floor
import numpy as np
import random 
import time
import itertools

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

def scaled_psi_4(s):
    '''
    Compute additional base theta function corresponding 
    to the case ZZ+1/4

    :param s: sampling width

    :returns: a real value
    '''
    q=nome(s)**4
    return (1/2)*jtheta(2, 0, q)

def theta_BW16(s):
    '''
    Compute the theta series of BW16

    :param s: sampling width

    :returns: a real value
    '''
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta4 = jtheta(4,0,q**4)
    theta_series = (1/2)*theta3**16+15*(theta3**8)*(theta2**8)+(1/2)*theta2**16+(1/2)*theta4**16
    return theta_series

def coset_prob_BW16(s):
    '''
    Compute the probability that a vector lies in a 
    coset of BW16 represented by a codeword of weight w=0,8,16

    :param s: sampling width

    :returns: a list of real values
    '''
    q = nome(s)
    theta3 = jtheta(3,0,q**4)
    theta2 = jtheta(2,0,q**4)
    theta4 = jtheta(4,0,q**4)
    theta_series = theta_BW16(s)
    table = []
    coset_0 = ((1/2)*theta3**16+(1/2)*theta4**16)/theta_series
    coset_8 = (15*(theta2**8)*theta3**8)/theta_series
    coset_16 = ((1/2)*theta2**16)/theta_series
    table.extend([coset_0, coset_8, coset_16])
    return table

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

def weight_enum_e(n, s):
    '''
    Compute even weight enumerator

    :param n: lattice dimension
    :param s: sampling width

    :returns: a real value
    '''
    sum = 0
    theta3 = jtheta(3,0,nome(s)**4)
    theta2 = jtheta(2,0,nome(s)**4)
    for i in range(floor((n-1)/2)):
        sum += comb(n-1, 2*i)*(theta3**(n-1-2*i))*(theta2**(2*i))
    return sum

def weight_enum_o(n, s):
    '''
    Compute odd weight enumerator

    :param n: lattice dimension
    :param s: sampling width

    :returns: a real value
    '''
    sum = 0
    theta3 = jtheta(3,0,nome(s)**4)
    theta2 = jtheta(2,0,nome(s)**4)
    for i in range(1,floor((n+1)/2)):
        sum += comb(n-1, 2*i-1)*(theta3**(n-2*i))*(theta2**(2*i-1))
    return sum

def sample_dn(n, s):
    ''''
    Sample a lattice vector from Dn 

    :param n: lattice dimension
    :param s: sampling width
    
    :returns: a lattice vector
    '''
    indices = [i for i in range(n)]
    m = random.choices(indices[1:n/2+1], weights=coset_prob_dn(n,s))
    subset = random.sample(indices, k=2*m[0])
    x = [None] * n
    for j in range(len(subset)):
        x[subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)()
    not_subset = list(set(indices)-set(subset))
    for j in range(len(not_subset)):
        x[not_subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))()
    return x

def sample_dn_special_shift(n, s, t):
    ''''
    Sample a lattice vector from Dn with a shift of the form (a,...,a) with a=1/2

    :param n: lattice dimension
    :param s: sampling width
    
    :returns: a lattice vector
    '''
    c = list(random.sample(list(codes.ParityCheckCode(GF(2), n-1)), k=1)[0])
    c_real = [float(i) for i in c]
    x = [None] * n
    for i in range(n):
        x[i] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=(c_real[i]+float(t[i]))/2)()
    return x

def sample_dn_bar(n, s):
    ''''
    Sample a lattice vector from Dn + (1,0,...,0) 

    :param n: lattice dimension
    :param s: sampling width
    
    :returns: a lattice vector
    '''
    indices = [i for i in range(n)]
    m = random.choices(indices[1:n/2+1], weights=coset_prob_dn(n,s))[0] # weight_ls is same as in normal Dn alg
    q = nome(s)
    if m < n/2:
        num_heads = (jtheta(2,0,q**4)**2)
        p_heads = num_heads/(num_heads+(2*m)/(n-2*m)*(jtheta(3,0,q**4)**2))
    else:
        p_heads = 0

    coin = random.choices(['heads', 'tails'], weights=[p_heads,1-p_heads])[0]

    x = [None] * n
    if coin == 'heads':
        subset = random.sample(indices[1:], k=2*m)+[0]
        for j in range(len(subset)):
            x[subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)()
        not_subset = list(set(indices)-set(subset))
        for j in range(len(not_subset)):
            x[not_subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))()
    else:
        subset = random.sample(indices[1:], k=2*m-1)
        for j in range(len(subset)):
            x[subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))),c=1/2)()
        not_subset = list(set(indices)-set(subset))
        for j in range(len(not_subset)):
            x[not_subset[j]] = 2 * DGaussZ(sigma=float(s)/(2*float(sqrt(2*pi))))()
    return x


# Precompute table of codewords of weight 8 in Reed-Muller code

C1 = codes.BinaryReedMullerCode(1,4)

weight_8 = []
for word in C1:
    if list(word).count(1) == 8:
        weight_8.append(list(word))

def sample_BW16(s, weights_prob):
    ''''
    Sample a lattice vector from BW16

    :param s: sampling width
    :param weights_prob: a list of the probabilities of sampling in each coset, computed offline by coset_prob_BW16
    
    :returns: a lattice vector
    '''
    w = random.choices([0,8,16], weights=weights_prob)[0]
    if w == 0: 
        samp = 2*np.array(sample_dn(16, s/4))
        return samp
    elif w == 16:
        samp = 2*np.array(sample_dn_special_shift(16, s/4, t=[1/2]*16))
        return samp
    else: 
        c = random.sample(weight_8, k=1)[0]

        # throw a biased coin with prob even, sample over even-even part, otherwise odd-odd part
        num_even = (2**7)*(scaled_psi_4(s)**8)*theta_dn(8, s)
        denom_even = num_even + (2**7)*(scaled_psi_4(s)**8)*(jtheta(2, 0, nome(s)**4)*weight_enum_e(8,s)+jtheta(3, 0, nome(s)**4)*weight_enum_o(8,s))
        p_even = num_even/denom_even

        decomp = random.choices(['even-even', 'odd-odd'], weights=[p_even,1-p_even])
        supp_J = [i for i in range(16) if c[i]==1] # set supports 
        not_supp = list(set([i for i in range(16)])-set(supp_J))

        if decomp == 'even-even': # sample from even-even part
            t1 = [1/2] * 8
            x = sample_dn_special_shift(8, s/4, t1) # double check if we need to scale this by 2?
        
            y = 2*np.array(sample_dn(8, s/4))
        else: 
            t1 = [3/2] + [1/2] * 7
            x = sample_dn_special_shift(8, s/4, t1)
        
            y = 2*np.array(sample_dn_bar(8, s/4))

        samp = [None] * 16
        for j in supp_J:
            samp[j] = x[supp_J.index(j)]
        for j in not_supp:
            samp[j] = y[not_supp.index(j)]

    return samp

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

def espitau_E8_shift(s,s_prime,shift): 
    ''''
    Sample a lattice vector from E8 with a shifted center

    :param s: smoothing parameter of D8
    :param s_prime: smoothing parameter of E8
    :param shift: shifted center
    
    :returns: a lattice vector
    '''
    b = random.sample([0,1], k=1)[0]
    bh = np.array([b*1/2 for i in range(8)])
    x = espitau_dn_shift(8, s,shift-bh)
    samp = bh+x
    prob = rho_vec(samp, s_prime)/rho_vec(samp, s)
    resp = random.choices(['accept', 'reject'], weights=[prob, 1-prob])[0]
    if resp == 'accept':
        return samp
    else:
        samp = espitau_E8_shift(s,s_prime,shift)
        return samp

# compute rotated basis of RBW8 used in doubling construction of BW16

R = np.array([[1,1,0,0,0,0,0,0],
              [1,-1,0,0,0,0,0,0],
              [0,0,1,1,0,0,0,0],
              [0,0,1,-1,0,0,0,0],
              [0,0,0,0,1,1,0,0],
              [0,0,0,0,1,-1,0,0],
              [0,0,0,0,0,0,1,1],
              [0,0,0,0,0,0,1,-1]])

B = np.array([[2, 0, 0, 0, 0, 0, 0, 0],
              [0, 2, 0, 0, 0, 0, 0, 0],
              [0, 0, 2, 0, 0, 0, 0, 0],
              [0, 0, 0, 2, 0, 0, 0, 0],
              [1, 1, 1, 0, 1, 0, 0, 0],
              [1, 0, 1, 1, 0, 1, 0, 0],
              [1, 0, 0, 1, 1, 0, 1, 0],
              [1, 0, 0, 0, 1, 1, 0, 1]])

RB=matrix(ZZ, np.dot(B,R).tolist()) # convert to sage matrix
RB=RB.echelon_form()

# compute the coset representatives for BW8/RBW8
A = []

for a1, a3, a4, a5, a6, a7 in itertools.product(range(2), range(2), range(2), range(4), range(2), range(4)):
    rep = a1*B[1,:]+a3*B[3,:]+a4*B[4,:]+a5*B[5,:]+a6*B[6,:]+a7*B[7,:]
    A.append(rep)

def espitau_BW16(s):
    ''''
    Sample a lattice vector from BW16 

    :param s: smoothing parameter of BW8
    
    :returns: a lattice vector
    '''
    alpha = np.array(random.sample(A, k=1)[0])
    m = np.append(alpha, alpha)
    u1 = np.dot(espitau_E8_shift(2.97, s, -m[0:8]), R) # sampling BW8 and then rotate to sample RBW8
    u2 = np.dot(espitau_E8_shift(2.97, s, -m[8:16]),R) 
    u = np.append(u1, u2)
    samp = [u[i] + m[i] for i in range(16)]
    return samp

'''
Compute the time it takes to sample 100000 times with both samplers and compare. 
Sampling width s is computed with epsilon = 2**(-36) and we take the example of n = 8
'''

# offline stage
coset_prob = coset_prob_BW16(2.3)

# online stage
start_time1 = time.perf_counter ()
for i in range(100000):
    sample_BW16(2.3, coset_prob)
end_time1 = time.perf_counter ()
time1 = end_time1 - start_time1

start_time2 = time.perf_counter ()
for i in range(100000):
    espitau_BW16(2.85)
end_time2 = time.perf_counter ()
time2 = end_time2 - start_time2

print("Our sampler is", time2/time1, "times as fast as the sampler by Espitau et. al.")

'''
Note: in order to obtain a valid probability prob = rho_vec(samp, s_prime)/rho_vec(samp, 2.97) we have that s_prime has to be at most 2.97.
The smoothing parameter of BW16 is about 2.3 for epsilon = 2^-36. However, when we tried to run the algorithm with 2.3 as input, we obtained error 
messages about reaching maximum depth recursion. We have to choose a sampling width above the smoothing of BW16, but smaller than 2.97 so we have a 
small available range.
'''