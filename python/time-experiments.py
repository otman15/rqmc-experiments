#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Timing experiment for the RQMC integrator using the built-in replication
parameter.
"""

import csv
import time
from pathlib import Path

import numpy as np
from qmcpy import DigitalNetB2

from integrands import mc2


def run_timing_experiment(k, s, m, seed, num_runs=3):
    n = 2**k
    timings = []

    for run in range(num_runs):
        start_time = time.perf_counter()

        estimates = np.empty(m)

        for replication in range(m):
            nus = DigitalNetB2(
                dimension=s,
                order="GRAY",
                randomize="NUS",
                replications=1,
                seed=seed + run + replication,
                t=32,
            )

            points = nus(n)
            estimates[replication] = np.mean(mc2(points))

        _ = estimates

        elapsed_time = time.perf_counter() - start_time
        timings.append(elapsed_time)

    return float(np.median(timings))


def main():
    s_values = [2, 4]
    k_values = [8, 10]
    m_values = [1, 5, 10, 50]
    seed = 20171215
    num_runs = 10

    output_dir = Path(__file__).resolve().parent / "py-results"
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    for s in s_values:
        for k in k_values:
            for m in m_values:
                median_time = run_timing_experiment(k, s, m, seed, num_runs)

                results.append(
                    {
                        "method": "QMCPy",
                        "s": s,
                        "k": k,
                        "m": m,
                        "median_time": median_time,
                    }
                )

                print(
                    f"s={s}, k={k}, m={m}, median time={median_time:.6f} seconds"
                )

    output_file = output_dir / "qmcpy_time_by_rep.csv"

    with output_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=["method", "s", "k", "m", "median_time"],
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"Results written to {output_file}")


if __name__ == "__main__":
    main()
