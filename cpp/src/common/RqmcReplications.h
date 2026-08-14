#ifndef RQMC_REPLICATIONS_H
#define RQMC_REPLICATIONS_H

#include <time.h>
#include <cstdint>
#include <string>
#include <vector>

struct RqmcReplicationResult
{
    std::vector<double> estimates;
    double elapsedTimeSeconds;
};

double evaluateRqmcModel(
    const std::string& model,
    const std::vector<double>& x);

template <typename Sequence>
RqmcReplicationResult runRqmcReplications(
    const std::string& model,
    unsigned int k,
    unsigned int s,
    unsigned int m,
    long seed)
{
    std::vector<double> estimates;
    estimates.reserve(m);

    std::vector<double> x(s);

    timespec start;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

    for (unsigned int r = 0; r < m; ++r)
    {
        Sequence sequence(
            s,
            k,
            seed + static_cast<long>(r));

        std::uint64_t n = sequence.numPoints();
        double sum = 0.0;

        for (std::uint64_t i = 0; i < n; ++i)
        {
            for (unsigned int j = 0; j < s; ++j)
                x[j] = sequence.coordinate(j);

            sum += evaluateRqmcModel(model, x);

            if (i + 1 < n)
                sequence.next();
        }

        estimates.push_back(
            sum / static_cast<double>(n));
    }

    timespec end;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    
    double elapsedTimeSeconds =
      static_cast<double>(end.tv_sec - start.tv_sec)
      + static_cast<double>(end.tv_nsec - start.tv_nsec) * 1e-9;

    return {estimates, elapsedTimeSeconds};
}

#endif
