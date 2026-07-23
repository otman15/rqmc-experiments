import csv
import time
from pathlib import Path

import numpy as np
from scipy.stats import kurtosis
from qmcpy import DigitalNetB2

from integrands import mc2


def run_rqmc(k, s, m, seed):
    n = 2**k
    estimates = np.empty(m)

    seeds = np.random.SeedSequence(seed).spawn(m)

    start_time = time.process_time()

    # With manual replications ##############################

    for randomization in range(m):
        nus = DigitalNetB2(
            dimension=s,
            randomize="NUS",
            seed=seeds[randomization]
        )

        points = nus(n)
        estimates[randomization] = np.mean(mc2(points))

    ########################################################

    '''
    ########################################################
    # Use built-in replications.
    # It is faster for small point sets, but for large ones
    # its execution time becomes similar to the manual loop.

    nus = DigitalNetB2(
        dimension=s,
        randomize="NUS",
        replications=m,
        seed=seed
    )

    points = nus(n)
    estimates = np.mean(mc2(points), axis=1)

    ########################################################
    '''

    cpu_time = time.process_time() - start_time

    '''
    mean_estimate = np.mean(estimates)
    centered = estimates - mean_estimate
    second_moment = np.mean(centered**2)
    kurtosis_numpy = (np.mean(centered**4) / second_moment**2 - 3.0)
    '''

    mean_estimate = np.mean(estimates)
    variance = np.var(estimates, ddof=1)
    kurtosis_scipy = kurtosis(estimates,fisher=True, bias=True)

    # print("With loop replications")
    print(f"dimension s              = {s}")
    print(f"k                        = {k}")
    print(f"points per randomization = {n}")
    print(f"randomizations m         = {m}")
    print(f"mean estimate            = {mean_estimate:.16e}")
    print(f"variance                 = {variance:.16e}")
    # print(f"kurtosis NumPy           = {kurtosis_numpy:.16e}")
    print(f"kurtosis SciPy           = {kurtosis_scipy:.16e}")
    print(f"CPU time                 = {cpu_time:.6f} seconds")
    print()

    return {
        "s": s,
        "k": k,
        "method": "QMCPy-NUS",
        "variance": variance,
        "cpu_time": cpu_time
    }


s_values = [2, 4]
k_values = [8 , 10, 12]

m = 10000
seed = 20171215

results = []

for s in s_values:
    for k in k_values:
        results.append(run_rqmc(k, s, m, seed))


output_file = Path("../results/qmcpy_results.csv")
output_file.parent.mkdir(parents=True, exist_ok=True)

with output_file.open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "s",
            "k",
            "method",
            "variance",
            "cpu_time"
        ]
    )

    writer.writeheader()
    writer.writerows(results)

print(f"Results written to {output_file}")