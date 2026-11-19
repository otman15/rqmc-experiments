#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "CachedNUSPointSet.h"
#include "FriedelKellerLibseq.h"
#include "HelmerStochasticSobol.h"
#include "NewCachedNusPointSet.h"
#include "PbrtSobolCore.h"
#include "RqmcReplications.h"
#include "Statistics.h"
#include "xor_values.h"

namespace
{
void verifyHelmerSobol(unsigned int s)
{
    if (s == 0 || s > static_cast<unsigned int>(sampling::MAX_SOBOL_DIM))
        throw std::invalid_argument("Common Sobol comparison requires 1..64 dimensions");

    // Helmer stores off-diagonal columns of the inverse Sobol matrix. Check
    // C * (e_c XOR xor[c]) = e_c for all 32 columns, outside measured work.
    for (unsigned int j = 0; j < s; ++j)
        for (unsigned int c = 0; c < 32; ++c)
            if (PbrtSobolCore::sobolInteger(
                    (std::uint64_t{1} << c) ^ sampling::sobol_xors[j][c], j)
                != (std::uint32_t{1} << (31 - c)))
                throw std::runtime_error("Helmer and PBRT Sobol matrices do not match");
}

std::vector<std::uint32_t> sobolBinaryColumns(unsigned int s, unsigned int k)
{
    std::vector<std::uint32_t> columns(static_cast<std::size_t>(s) * k);
    for (unsigned int j = 0; j < s; ++j)
        for (unsigned int c = 0; c < k; ++c)
            columns[static_cast<std::size_t>(j) * k + c] =
                PbrtSobolCore::sobolInteger(std::uint64_t{1} << c, j);
    return columns;
}

template <typename DigitalNet>
DigitalNet commonSobolNet(unsigned int s, unsigned int k)
{
    const auto binaryColumns = sobolBinaryColumns(s, k);
    DigitalNet net{
        std::uint64_t{1} << k, s, k, 32,
        std::vector<std::uint64_t>(static_cast<std::size_t>(s) * k)};

    for (unsigned int j = 0; j < s; ++j)
    {
        std::uint64_t grayColumn = 0;
        for (unsigned int c = 0; c < k; ++c)
        {
            // Both cached implementations enumerate Gray-code indices. The
            // cumulative XOR converts their basis to the binary-index order
            // used by Helmer's tables and libseq: D * Gray(i) = C * i.
            grayColumn ^= binaryColumns[static_cast<std::size_t>(j) * k + c];
            net.columns[static_cast<std::size_t>(j) * k + c] = grayColumn;
        }
    }
    return net;
}

struct MethodResult
{
    std::string name;
    std::vector<double> estimates;
    double setupSeconds = 0.0;
    double scramblingSeconds = 0.0;
    double integrationSeconds = 0.0;
    bool constructorScrambles = false;
};

template <typename Sequence>
double integrate(Sequence& sequence, const std::string& model, std::vector<double>& x)
{
    sequence.reset();
    const auto n = sequence.numPoints();
    double sum = 0.0;
    for (std::uint64_t i = 0; i < n; ++i)
    {
        for (unsigned int j = 0; j < x.size(); ++j)
            x[j] = sequence.coordinate(j);
        sum += evaluateRqmcModel(model, x);
        if (i + 1 < n)
            sequence.next();
    }
    return sum / static_cast<double>(n);
}

template <typename PrepareSequence>
MethodResult runMethod(
    const std::string& name,
    PrepareSequence prepare,
    const std::string& model,
    unsigned int s,
    unsigned int m,
    long seed,
    unsigned int warmup,
    bool constructorScrambles = false)
{
    MethodResult result;
    result.name = name;
    result.constructorScrambles = constructorScrambles;
    // Reserve common bookkeeping buffers outside the measured passes.
    result.estimates.reserve(m);
    std::vector<double> x(s);

    // Measure the constructor's complete setup work. Both NUS variants and HCK
    // prepare only; libseq's mandatory initial scrambling is included and reported.
    auto start = RqmcDetail::processCpuTime();
    auto sequence = prepare();
    result.setupSeconds = RqmcDetail::secondsBetween(
        start, RqmcDetail::processCpuTime());

    for (unsigned int r = 0; r < warmup; ++r)
    {
        sequence.randomize();
        (void) integrate(sequence, model, x);
    }

    // Seed once per measured pass and advance the same RNG across replications.
    // Reseeding discards warmup's RNG advancement without rebuilding anything.
    sequence.reseed(static_cast<std::uint64_t>(seed));
    start = RqmcDetail::processCpuTime();
    for (unsigned int r = 0; r < m; ++r)
    {
        sequence.randomize();
        result.estimates.push_back(integrate(sequence, model, x));
    }
    result.integrationSeconds = RqmcDetail::secondsBetween(
        start, RqmcDetail::processCpuTime());

    // A separate pass, replaying the same stream, measures cache generation
    // without function evaluation. It reuses the same prepared object/storage.
    sequence.reseed(static_cast<std::uint64_t>(seed));
    start = RqmcDetail::processCpuTime();
    for (unsigned int r = 0; r < m; ++r)
        sequence.randomize();
    result.scramblingSeconds = RqmcDetail::secondsBetween(
        start, RqmcDetail::processCpuTime());
    return result;
}

void printReport(const std::vector<MethodResult>& results, unsigned int m)
{
    std::cout << "\nStatistics of the " << m << " measured integral estimates.\n"
              << "Skewness and excess kurtosis are bias corrected.\n"
              << std::left << std::setw(26) << "Method" << std::right
              << std::setw(16) << "Mean"
              << std::setw(16) << "Variance"
              << std::setw(16) << "Skewness"
              << std::setw(20) << "Excess kurtosis" << '\n';
    std::cout << std::scientific << std::setprecision(6);
    for (const auto& result : results)
        std::cout << std::left << std::setw(26) << result.name << std::right
                  << std::setw(16) << mean(result.estimates)
                  << std::setw(16) << sampleVariance(result.estimates)
                  << std::setw(16) << skewness(result.estimates)
                  << std::setw(20) << excessKurtosis(result.estimates) << '\n';

    std::cout << "\nProcess CPU seconds; scrambling and integration are totals over " << m << " replications.\n"
              << "Setup includes construction, allocations, preprocessing, and any automatic initial scrambling.\n"
              << "Total = setup + scrambling. Integration includes scrambling and function evaluation.\n"
              << "Integration excludes setup and warmup. Both measured passes start from the printed seed.\n"
              << std::left << std::setw(26) << "Method" << std::right
              << std::setw(16) << "Setup"
              << std::setw(16) << "Scrambling"
              << std::setw(16) << "Total"
              << std::setw(20) << "Integration" << '\n';
    std::cout << std::fixed << std::setprecision(6);
    for (const auto& result : results)
        std::cout << std::left << std::setw(26) << result.name << std::right
                  << std::setw(16) << result.setupSeconds
                  << std::setw(16) << result.scramblingSeconds
                  << std::setw(16) << result.setupSeconds + result.scramblingSeconds
                  << std::setw(20) << result.integrationSeconds << '\n';
    for (const auto& result : results)
        if (result.constructorScrambles)
            std::cout << result.name << ": setup includes upstream's mandatory initial scrambling; "
                      << "each measured pass still performs " << m << " fresh randomizations.\n";
}
}

int main()
{
    // Edit these settings, rebuild compare_nus_methods, then run it again.
    // Models: gaussian, mc2, sumueu, smoothperb4, polynomial, oscillatory,
    //         smoothgauss, piecelingauss, indsumnormal.
    const std::string model = "mc2";
    const unsigned int s = 10;        // Dimension (HCK supports 1..64).
    const unsigned int k = 18;        // Number of points n = 2^k; k < 31.
    const unsigned int m = 1000;       // Measured replications.
    const long seed = 20171215;
    const unsigned int warmup = 100;   // Discarded replications; zero disables warmup.
    // Variance requires m >= 2, skewness m >= 3, and excess kurtosis m >= 4;
    // the existing statistics functions print NaN when these are undefined.

    try
    {
        if (s == 0 || k >= 31 || m == 0)
            throw std::invalid_argument("Require s > 0, k < 31, and m > 0");
        if (seed < 0 || static_cast<std::uint64_t>(seed) >
                            std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("Common RNG comparison requires a seed in 0..2^32-1");
        verifyHelmerSobol(s);
        std::cout << "model=" << model << ", s=" << s << ", k=" << k
                  << ", n=" << (std::uint64_t{1} << k) << ", m=" << m
                  << ", warmup=" << warmup << ", seed=" << seed << '\n';
        std::cout << "Experimental common inputs: Helmer PCG32 and PBRT Sobol matrices "
                     "in binary-index order, verified against Helmer's tables.\n"
                     "Libseq and CachedNUSPointSet receive PCG32 through their existing drand48 interface.\n"
                     "Helmer retains native within-stratum jitter; the other methods use 32 output bits.\n";

        std::vector<MethodResult> results;
        results.reserve(4);
        // Comment/uncomment a method call to select which methods to run.


        results.push_back(runMethod("HelmerStochasticSobol", [=]
        {
            return HelmerStochasticSobol(s, k, seed,
                HelmerStochasticSobol::Initialization::PrepareOnly);
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("NewCachedNusPointSet", [=]
        {
            return NewCachedNusPointSet(
                commonSobolNet<NewCachedNusPointSet::DigitalNet>(s, k),
                static_cast<std::uint64_t>(seed));
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("CachedNUSPointSet", [=]
        {
            return CachedNUSPointSet(commonSobolNet<CachedNUSPointSet::DigitalNet>(s, k),
                static_cast<std::uint64_t>(seed));
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("Friedel-Keller/libseq", [=]
        {
            return FriedelKellerLibseq(s, k, seed, sobolBinaryColumns(s, k));
        }, model, s, m, seed, warmup, true));

        printReport(results, m);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Comparison failed: " << error.what() << '\n';
        return 1;
    }
}
