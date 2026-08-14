from pathlib import Path
import time

import numpy as np
from qmcpy import DigitalNetB2

from integrands import mc2

'''
the goal of this is to test if using the persistent qmcpy nus can save time in simulations.
by persistent we mean using the alredy generated points (2*k) when to get 2*(k+i) points
'''
def run_persistent_rqmc(
    model_name, method_name,
    s_values, k_values, m, seed,
    output_dir,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ordered_k_values = sorted(k_values)

    for s in s_values:
        estimates_by_k = {
            k: np.empty(m)
            for k in ordered_k_values
        }
        incremental_times_by_k = {
            k: 0.0
            for k in ordered_k_values
        }

        for replication in range(m):
            replication_start_time = time.perf_counter()

            # One fixed NUS mapping for this replication, reused at every k.
            sobol = DigitalNetB2(
                dimension=s,
                randomize="NUS",
                order="GRAY",
                replications=1,
                seed=seed + replication,
                t=32,
            )

            cumulative_sum = 0.0
            previous_n = 0

            for k_index, k in enumerate(ordered_k_values):
                # Setup is part of the initial estimate only. Later estimates
                # time just the newly added point block and its evaluation.
                increment_start_time = (
                    replication_start_time
                    if k_index == 0
                    else time.perf_counter()
                )
                current_n = 1 << k

                # Generate only the points added for this k.
                new_points = sobol(
                    n_min=previous_n,
                    n_max=current_n,
                )

                # With one replication, the result has shape (1, new_n).
                cumulative_sum += np.sum(mc2(new_points))
                estimates_by_k[k][replication] = cumulative_sum / current_n
                incremental_times_by_k[k] += (
                    time.perf_counter() - increment_start_time
                )

                previous_n = current_n

        for k in ordered_k_values:
            filename = (
                output_dir
                / f"{model_name}-{s}-{method_name}-{k}-{m}.dat"
            )

            # One sorted estimate per line, as in the other experiments.
            np.savetxt(filename, np.sort(estimates_by_k[k]), fmt="%.17g")
            print(
                f"model={model_name}, s={s}, k={k}, m={m}, "
                f"incremental_time={incremental_times_by_k[k]:.6f} seconds"
            )
            print(f"Results written to {filename}")


if __name__ == "__main__":
    experiment_start_time = time.perf_counter()

    model_name = "MC2"
    method_name = "NUS-QMCPy"

    s_values = [4]
    k_values = [8, 10, 12, 14]
    m = 1000
    seed = 20171215

    output_dir = Path(__file__).resolve().parent / "py-results" / "dat2"

    run_persistent_rqmc(
        model_name,
        method_name,
        s_values,
        k_values,
        m,
        seed,
        output_dir,
    )
    print(
        "Total experiment time="
        f"{time.perf_counter() - experiment_start_time:.6f} seconds"
    )
