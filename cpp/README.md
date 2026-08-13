# C++ NUS Experiments

C++ implementations used to compare Nested Uniform Scrambling / Owen scrambling methods.

# C++ NUS Experiments

C++ implementations used to compare Nested Uniform Scrambling / Owen scrambling methods.

## Dependencies

Required system packages:

- CMake
- C++17-compatible compiler
- Boost development libraries

Ubuntu/Debian:

```bash
sudo apt update
sudo apt install build-essential cmake libboost-all-dev
```
## Implementations

### Friedel–Keller / libseq

Original C++ implementation by Friedel and Keller for randomized low-discrepancy point sets.

Source:  
https://multires.caltech.edu/software/libseq/

### Burley / QuantLib

QuantLib implementation of Brent Burley's Practical Hash-based Owen Scrambling.

Source:  
https://github.com/lballabio/QuantLib

Paper:  
https://jcgt.org/published/0009/04/01/

### PBRT OwenScrambler

PBRT-v4 implementation of prefix-dependent hashed Owen scrambling.

Source:  
https://github.com/mmp/pbrt-v4

### PBRT FastOwenScrambler

PBRT-v4 fast approximate Owen scrambling implementation.

Source:  
https://github.com/mmp/pbrt-v4

### Helmer–Christensen–Kensler Stochastic Sobol

C++ implementation accompanying *Stochastic Generation of (t,s) Sample Sequences*.

Source:  
https://github.com/Andrew-Helmer/stochastic-generation

### Adaptive Regular Tiles / ART-Owen

Official implementation of ART-Owen scrambling based on Adaptive Regular Tiles.

Source:  
https://abdallagafar.com/publications/artowen/

Its integration into the common C++ experiments is currently unfinished.

## Build

From the `cpp/` directory:

```bash
cmake -S . -B build
cmake --build build -j$(nproc)
```

## Timing Experiment

Run:

```bash
./build/time_experiments
```

The experiment runs all currently integrated C++ methods.

Timing parameters are defined in:

`src/experiments/TimeExperiments.cpp`

In `main()`, edit:

- `models`
- `sValues`
- `kValues`
- `timingMValues`
- `seed`
- `numRuns`

Timing results are written to:

`cpp-results/`

## Statistics Experiment

Run:

```bash
./build/statistics_experiments
```

The experiment runs all currently integrated C++ methods.

Statistics parameters are defined in:

`src/experiments/StatisticsExperiments.cpp`

In `main()`, edit:

- `models`
- `sValues`
- `kValues`
- `statisticsM`
- `seed`

Statistics samples are written to:

`cpp-results/dat/`
