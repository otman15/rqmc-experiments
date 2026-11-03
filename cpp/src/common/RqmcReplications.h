#ifndef RQMC_REPLICATIONS_H
#define RQMC_REPLICATIONS_H

#include <time.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

struct RqmcReplicationResult
{
    std::vector<double> estimates;
    double elapsedTimeSeconds;
    // Populated by the prepared runner only. All values are process CPU seconds.
    double preparationTimeSeconds = 0.0;
    double scramblingTimeSeconds = 0.0;
    double evaluationTimeSeconds = 0.0;
};

double evaluateRqmcModel(
    const std::string& model,
    const std::vector<double>& x);

namespace RqmcDetail
{
inline timespec processCpuTime()
{
    timespec time;
    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time) != 0)
        throw std::runtime_error("Cannot read the process CPU clock");
    return time;
}

inline double secondsBetween(const timespec& start, const timespec& end)
{
    return static_cast<double>(end.tv_sec - start.tv_sec) +
           static_cast<double>(end.tv_nsec - start.tv_nsec) * 1e-9;
}
}

// prepare() must construct and seed a sequence without scrambling it. Unlike
// the legacy runner below (seed + replication, new sequence each time), this
// path seeds once per run and advances the same RNG across all replications.
// The factory also supports nonmovable prepared sequences through C++17 elision.
//
// Preparation includes the estimate/model buffers and all sequence allocations
// and preprocessing. Scrambling totals include randomize() and the rewind.
// Evaluation totals include reading coordinates, iteration, model calls, the
// sample mean, and storing the estimate. elapsedTimeSeconds covers the complete
// experiment, including preparation, timer/bookkeeping overhead and cleanup;
// it excludes reporting/statistics performed by the caller. Phase totals are
// over all m replications, not per-replication averages.
template <typename PrepareSequence>
RqmcReplicationResult runPreparedRqmcReplications(
    const std::string& model,
    unsigned int s,
    unsigned int m,
    PrepareSequence prepare)
{
    RqmcReplicationResult result;
    const auto experimentStart = RqmcDetail::processCpuTime();
    {
        result.estimates.reserve(m);
        std::vector<double> x(s);
        auto sequence = prepare();
        const std::uint64_t n = sequence.numPoints();
        if (n == 0)
            throw std::invalid_argument("RQMC sequence must contain points");
        result.preparationTimeSeconds = RqmcDetail::secondsBetween(
            experimentStart, RqmcDetail::processCpuTime());

        for (unsigned int r = 0; r < m; ++r)
        {
            const auto scrambleStart = RqmcDetail::processCpuTime();
            sequence.randomize();
            sequence.reset();
            const auto evaluationStart = RqmcDetail::processCpuTime();
            result.scramblingTimeSeconds += RqmcDetail::secondsBetween(
                scrambleStart, evaluationStart);

            double sum = 0.0;
            for (std::uint64_t i = 0; i < n; ++i)
            {
                for (unsigned int j = 0; j < s; ++j)
                    x[j] = sequence.coordinate(j);
                sum += evaluateRqmcModel(model, x);
                if (i + 1 < n)
                    sequence.next();
            }
            result.estimates.push_back(sum / static_cast<double>(n));
            result.evaluationTimeSeconds += RqmcDetail::secondsBetween(
                evaluationStart, RqmcDetail::processCpuTime());
        }
    }
    result.elapsedTimeSeconds = RqmcDetail::secondsBetween(
        experimentStart, RqmcDetail::processCpuTime());
    return result;
}

// Legacy lifecycle retained for adapters that generate points on construction.
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
