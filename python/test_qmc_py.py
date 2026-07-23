import time
import numpy as np
from scipy.stats import kurtosis
from qmcpy import DigitalNetB2

from integrands import mc2


def run_rqmc(k, s, m, seed):
    n = 2**k
    estimates = np.empty(m)

    seeds = np.random.SeedSequence(seed).spawn(m)

    start_time = time.perf_counter()
    
    #with manual replications##############################
    
    for randomization in range(m):
        nus = DigitalNetB2(
            dimension=s,
            randomize="NUS",
            seed=seeds[randomization]
        )

        points = nus(n)
        estimates[randomization] = np.mean(mc2(points))
    #################################
    '''
    # use built-in replications######################
    nus = DigitalNetB2(
        dimension=s,
        randomize="NUS",
        replications=m,
        seed=seed
    )
    
    points = nus(n)
    estimates = np.mean(mc2(points), axis=1)
    ##########################################  
    '''
    
    elapsed_time = time.perf_counter() - start_time
      
    mean_estimate = np.mean(estimates)
    centered = estimates - mean_estimate
    second_moment = np.mean(centered**2)
    
    kurtosis_numpy = np.mean(centered**4) / second_moment**2
    kurtosis_scipy = kurtosis(estimates, fisher=False, bias=True)
    
    #print("With loop replications")
    print(f"dimension s              = {s}")
    print(f"k                        = {k}")
    print(f"points per randomization = {n}")
    print(f"randomizations m         = {m}")
    print(f"mean estimate            = {np.mean(estimates):.16e}")
    print(f"variance                 = {np.var(estimates, ddof=1):.16e}")
    print(f"kurtosis NumPy          = {kurtosis_numpy:.16e}")
    print(f"kurtosis SciPy          = {kurtosis_scipy:.16e}")
    print(f"total time               = {elapsed_time:.6f} seconds")
    print()

    return estimates


s_values = [2, 4]
k_values = [4, 8, 10, 12]

m = 10000
seed = 20171215

for s in s_values:
    for k in k_values:
        run_rqmc(k, s, m, seed)