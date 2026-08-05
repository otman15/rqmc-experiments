#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RQMC timing or statistics experiment using manual independent replications.
"""

import csv
import time
from functools import partial
from pathlib import Path

import numpy as np
from qmcpy import DigitalNetB2
from scipy.stats import kurtosis, skew

from integrands import mc2, polynomial, smoothperb4, sumueu


def run_rqmc_replications(model, k, s, m, seed):
    """Run m RQMC replications once and return their estimates and elapsed time."""
    n = 2**k
    start_time = time.perf_counter()
    estimates = np.empty(m)

    for replication in range(m):
        nus = DigitalNetB2(
            dimension=s,
            order="GRAY",
            randomize="NUS",
            replications=1,
            seed=seed + replication,
            t=32,
        )

        points = nus(n)
        estimates[replication] = np.mean(model(points))

    elapsed_time = time.perf_counter() - start_time
    return estimates, elapsed_time


def run_timing_experiments(
    models, s_values, k_values, timing_m_values, seed, num_runs, output_dir
):
    output_dir.mkdir(parents=True, exist_ok=True)

    for model_name, model in models:
        results = []

        for s in s_values:
            for k in k_values:
                for m in timing_m_values:
                    timings = []

                    for run in range(num_runs):
                        _, elapsed_time = run_rqmc_replications(
                            model, k, s, m, seed + run*m
                        )
                        timings.append(elapsed_time)

                    median_time = float(np.median(timings))
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
                        f"model={model_name}, s={s}, k={k}, m={m}, "
                        f"median time={median_time:.6f} seconds"
                    )

        output_file = output_dir / f"qmcpy_{model_name}_time_rep-{num_runs}.csv"

        with output_file.open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(
                file,
                fieldnames=["method", "s", "k", "m", "median_time"],
            )
            writer.writeheader()
            writer.writerows(results)

        print(f"Results written to {output_file}")


def run_statistics_experiments(
    models, s_values, k_values, statistics_m, seed, output_dir, method
):
    output_dir.mkdir(parents=True, exist_ok=True)

    for model_name, model in models:
        for s in s_values:
            for k in k_values:
                estimates, elapsed_time = run_rqmc_replications(
                    model, k, s, statistics_m, seed
                )
                variance = (
                    np.var(estimates, ddof=1) if statistics_m > 1 else np.nan
                )
                skewness = (
                    skew(estimates, bias=False) if statistics_m > 2 else np.nan
                )
                excess_kurtosis = (
                    kurtosis(estimates, fisher=True, bias=False)
                    if statistics_m > 3
                    else np.nan
                )

                print(
                    f"model={model_name}, s={s}, k={k}, m={statistics_m}, "
                    f"mean={np.mean(estimates):.16e}, "
                    f"variance={variance:.16e}, "
                    f"skewness={skewness:.16e}, "
                    f"excess kurtosis={excess_kurtosis:.16e}, "
                    f"time={elapsed_time:.6f} seconds"
                )

                output_file = output_dir / (
                    f"{model_name}-{s}-{method}-{k}-{statistics_m}.dat"
                )
                np.savetxt(output_file, np.sort(estimates), fmt="%.17g")
                print(f"Results written to {output_file}")


def main():
    timing_output_dir = Path(__file__).resolve().parent / "py-results" #path("path/")
    statistics_output_dir = timing_output_dir / "dat"
    TIMING = "TIMING"
    STATISTICS = "STATISTICS"
    mode = STATISTICS

    method = "NUS-QMCPy"
    models = [
        ("MC2", mc2),
        ("SumUeU", sumueu),
        ("SmoothPerB4", partial(smoothperb4, expw=1.0)),
        ("Polynomial", polynomial),
    ]
    s_values = [2, 4, 8, 16, 32]
    k_values = [8, 10, 12, 14, 16]
    timing_m_values = [1, 2, 5, 10, 20, 30, 50]
    statistics_m = 1000
    seed = 20171215
    num_runs = 11

    if mode == TIMING:
        run_timing_experiments(
            models,
            s_values,
            k_values,
            timing_m_values,
            seed,
            num_runs,
            timing_output_dir,
        )
    elif mode == STATISTICS:
        run_statistics_experiments(
            models,
            s_values,
            k_values,
            statistics_m,
            seed,
            statistics_output_dir,
            method,
        )
    else:
        raise ValueError(f"mode must be {TIMING!r} or {STATISTICS!r}, not {mode!r}")


if __name__ == "__main__":
    main()
