#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "C_matrix.h"
#include "CachedNUSPointSet.h"
#include "FriedelKellerLibseq.h"
#include "HelmerNativeSobol.h"
#include "HelmerStochasticSobol.h"
#include "NewCachedNusPointSet.h"
#include "RqmcReplications.h"
#include "Statistics.h"

namespace
{
CachedNUSPointSet::DigitalNet libseqNiederreiterNet(unsigned int s, unsigned int k)
{
    if (s == 0 || k >= 31)
        throw std::invalid_argument("Require s > 0 and k < 31 for the Niederreiter net");

    Ring ring(2);
    char generator[] = "Niederreiter";
    C_matrix matrix(&ring, generator, s);
    constexpr unsigned int outputBits = 32;
    CachedNUSPointSet::DigitalNet net{
        std::uint64_t{1} << k, s, k, outputBits,
        std::vector<std::uint64_t>(static_cast<std::size_t>(s) * k)};

    for (unsigned int j = 0; j < s; ++j)
    {
        std::uint64_t grayColumn = 0;
        for (unsigned int c = 0; c < k; ++c)
        {
            // Pack rows from most to least significant bit, exactly as
            // libseq's DigitalSequence_base_2::transfer_C() does.
            std::uint64_t binaryColumn = 0;
            for (unsigned int row = 0; row < outputBits; ++row)
            {
                binaryColumn <<= 1;
                if (row < matrix.NumDigits && c < matrix.NumDigits)
                    binaryColumn |= matrix.query(j, row, c);
            }

            // CachedNUSPointSet enumerates Gray-code indices. The XOR of
            // binary columns 0..c converts its basis to libseq's binary-index
            // point order: D * Gray(i) = C * i.
            grayColumn ^= binaryColumn;
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
    bool nativeSeeding = false;
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
    constexpr bool nativeSeeding = std::is_same_v<decltype(sequence), HelmerNativeSobol>;
    result.nativeSeeding = nativeSeeding;

    for (unsigned int r = 0; r < warmup; ++r)
    {
        sequence.randomize();
        (void) integrate(sequence, model, x);
    }

    // Seeded adapters replay the printed seed; the native function initializes
    // its own RNG on every call and exposes no seeding API.
    if constexpr (!nativeSeeding)
        sequence.reseed(static_cast<std::uint64_t>(seed));
    start = RqmcDetail::processCpuTime();
    for (unsigned int r = 0; r < m; ++r)
    {
        sequence.randomize();
        result.estimates.push_back(integrate(sequence, model, x));
    }
    result.integrationSeconds = RqmcDetail::secondsBetween(
        start, RqmcDetail::processCpuTime());

    // A separate pass measures generation without function evaluation. Seeded
    // adapters replay the integration stream; native calls use fresh seeds.
    if constexpr (!nativeSeeding)
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
              << "Integration excludes setup and warmup. Seeded methods replay the printed seed in both passes.\n"
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
    {
        if (result.constructorScrambles)
            std::cout << result.name << ": setup includes upstream's mandatory initial scrambling; "
                      << "each measured pass still performs " << m << " fresh randomizations.\n"
                      << "This method uses libseq's Niederreiter net and native drand48 RNG.\n";
        if (result.nativeSeeding)
            std::cout << result.name << ": directly calls upstream GetStochasticSobolSamples "
                      << "with shuffle=false, candidates=1, and owen=false.\n"
                      << "Its native PCG32 seeds itself on every call; the printed seed does not apply.\n"
                      << "Scrambling includes native RNG initialization and temporary allocations; "
                      << "the output buffer is reused. The two measured passes use different random samples.\n";
    }
}
}

int main()
{
    // Edit these settings, rebuild compare_nus_methods, then run it again.
    // Models: gaussian, mc2, sumueu, smoothperb4, polynomial, oscillatory,
    //         smoothgauss, piecelingauss, indsumnormal.
    const std::string model = "mc2";
    const unsigned int s = 10;        // Dimension (HCK supports 1..64).
    const unsigned int k = 20;        // Number of points n = 2^k; k < 31.
    const unsigned int m = 100;       // Measured replications.
    const long seed = 20171215;
    const unsigned int warmup = 100;   // Discarded replications; zero disables warmup.
    // Variance requires m >= 2, skewness m >= 3, and excess kurtosis m >= 4;
    // the existing statistics functions print NaN when these are undefined.

    try
    {
        if (s == 0 || k >= 31 || m == 0)
            throw std::invalid_argument("Require s > 0, k < 31, and m > 0");
        std::cout << "model=" << model << ", s=" << s << ", k=" << k
                  << ", n=" << (std::uint64_t{1} << k) << ", m=" << m
                  << ", warmup=" << warmup << ", seed=" << seed << '\n';
        std::cout << "CachedNUSPointSet uses libseq's 32-bit Niederreiter net "
                     "with matching point order.\n";

        std::vector<MethodResult> results;
        results.reserve(5);
        // Comment/uncomment a method call to select which methods to run.


        results.push_back(runMethod("HelmerStochasticSobol", [=]
        {
            return HelmerStochasticSobol(s, k, seed,
                HelmerStochasticSobol::Initialization::PrepareOnly);
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("Helmer-native", [=]
        {
            return HelmerNativeSobol(s, k);
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("NewCachedNusPointSet", [=]
        {
            return NewCachedNusPointSet(s, k, seed); // Prepares only; 32 output bits.
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("CachedNUSPointSet", [=]
        {
            return CachedNUSPointSet(libseqNiederreiterNet(s, k),
                static_cast<std::uint64_t>(seed));
        }, model, s, m, seed, warmup));

        results.push_back(runMethod("Friedel-Keller/libseq", [=]
        {
            return FriedelKellerLibseq(s, k, seed);
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
