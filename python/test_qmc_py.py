#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 12:16:30 2026

@author: otman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 15:17:58 2026

@author: otman
"""
import time
import numpy as np
from qmcpy import DigitalNetB2


def integrand(points):
    return np.prod(points, axis=1)


def run_rqmc(k, s, m, seed):
    n = 2**k
    estimates = np.empty(m)

    seeds = np.random.SeedSequence(seed).spawn(m)

    start_time = time.perf_counter()

    for randomization in range(m):
        nus = DigitalNetB2(
            dimension=s,
            randomize="NUS",
            seed=seeds[randomization]
        )

        points = nus(n)
        estimates[randomization] = np.mean(integrand(points))

    elapsed_time = time.perf_counter() - start_time

    print(f"With loop replications")
    print(f"dimension s              = {s}")
    print(f"k                        = {k}")
    print(f"points per randomization = {n}")
    print(f"randomizations m         = {m}")
    print(f"mean estimate            = {np.mean(estimates):.16e}")
    print(f"variance                 = {np.var(estimates, ddof=1):.16e}")
    print(f"total time               = {elapsed_time:.6f} seconds")

    return estimates

s = 4
k = 8
m = 100
seed = 20171215

estimates = run_rqmc(k, s, m, seed)
