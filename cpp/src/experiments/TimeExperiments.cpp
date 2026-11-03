#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "BurleyQuantLib.h"
#include "FriedelKellerLibseq.h"
#include "HelmerStochasticSobol.h"
#include "NewCachedNusPointSet.h"
#include "PbrtFastOwenScrambler.h"
#include "PbrtOwenScrambler.h"
#include "AdaptiveRegularTiles.h"
#include "RqmcReplications.h"
#include "Statistics.h"

enum class Method
{
    FriedelKellerLibseq,
    BurleyQuantLib,
    PbrtOwen,
    PbrtFastOwen,
    HelmerStochasticSobol
};

template <typename Sequence, bool ReuseStorage = false>
void runTimingForMethod(
    const std::string& methodName,
    const std::string& fileTag,
    const std::vector<std::pair<std::string, std::string>>& models,
    const std::vector<unsigned int>& sValues,
    const std::vector<unsigned int>& kValues,
    const std::vector<unsigned int>& timingMValues,
    long seed,
    unsigned int numRuns,
    const std::filesystem::path& outputDir)
{
    if constexpr (ReuseStorage)
        std::cout << methodName
                  << ": seed once per run; continuous RNG across replications. "
                  << "CPU phase totals cover all m replications; medians are over "
                  << numRuns << " runs and are computed independently.\n";
    for (const auto& [modelName, model] : models)
    {
        const std::filesystem::path outputFile =
            outputDir /
            (fileTag + "_" +
             modelName +
             "_time_rep-" +
             std::to_string(numRuns) +
             ".csv");

        std::ofstream file(outputFile);

        if (!file)
            throw std::runtime_error(
                "Cannot open output file: " +
                outputFile.string());

        if constexpr (ReuseStorage)
            file << "method,s,k,m,median_preparation_cpu_s,median_scrambling_total_cpu_s,"
                    "median_evaluation_total_cpu_s,median_experiment_total_cpu_s\n";
        else
            file << "method,s,k,m,median_time\n";
        file << std::setprecision(17);

        for (unsigned int s : sValues)
        {
            for (unsigned int k : kValues)
            {
                for (unsigned int m : timingMValues)
                {
                    std::vector<double> timings;
                    timings.reserve(numRuns);
                    std::vector<double> preparationTimes, scramblingTimes, evaluationTimes;
                    if constexpr (ReuseStorage)
                    {
                        preparationTimes.reserve(numRuns);
                        scramblingTimes.reserve(numRuns);
                        evaluationTimes.reserve(numRuns);
                    }

                    for (unsigned int run = 0;
                         run < numRuns;
                         ++run)
                    {
                        const long runSeed =
                            seed +
                            static_cast<long>(run) *
                            static_cast<long>(m);

                        const auto result = [&]
                        {
                            if constexpr (ReuseStorage)
                                return runPreparedRqmcReplications(model, s, m, [=]
                                {
                                    if constexpr (std::is_same_v<Sequence, HelmerStochasticSobol>)
                                        return Sequence(s, k, runSeed, Sequence::Initialization::PrepareOnly);
                                    else
                                        return Sequence(s, k, runSeed);
                                });
                            else
                                return runRqmcReplications<Sequence>(model, k, s, m, runSeed);
                        }();

                        timings.push_back(
                            result.elapsedTimeSeconds);
                        if constexpr (ReuseStorage)
                        {
                            preparationTimes.push_back(result.preparationTimeSeconds);
                            scramblingTimes.push_back(result.scramblingTimeSeconds);
                            evaluationTimes.push_back(result.evaluationTimeSeconds);
                        }
                    }

                    const double medianTime =
                        median(timings);

                    file
                        << methodName << ','
                        << s << ','
                        << k << ','
                        << m << ',';

                    if constexpr (ReuseStorage)
                        file << median(preparationTimes) << ','
                             << median(scramblingTimes) << ','
                             << median(evaluationTimes) << ',';
                    file << medianTime << '\n';

                    std::cout
                        << "method=" << methodName
                        << ", model=" << modelName
                        << ", s=" << s
                        << ", k=" << k
                        << ", m=" << m
                        << std::fixed
                        << std::setprecision(6);
                    if constexpr (ReuseStorage)
                        std::cout << ", median preparation CPU(s)=" << median(preparationTimes)
                                  << ", median scrambling CPU(s), total over m=" << median(scramblingTimes)
                                  << ", median evaluation CPU(s), total over m=" << median(evaluationTimes)
                                  << ", median experiment CPU(s), total=" << medianTime << '\n';
                    else
                        std::cout << ", median time=" << medianTime << " seconds\n";
                }
            }
        }

        std::cout
            << "Results written to "
            << outputFile << '\n';
    }
}

int main()
{
    const std::filesystem::path outputDir = "cpp-results";
    // const std::filesystem::path outputDir = "/mnt/c/Users/cherrato/Dropbox/Nus-comparisons/all-methods";
    std::filesystem::create_directories(outputDir);

  const std::vector<std::pair<std::string, std::string>> models = {
      {"MC2", "mc2"}
  };

  const std::vector<unsigned int> sValues = {
      2, 4, 8, 16
  };

  const std::vector<unsigned int> kValues = {
      8, 10, 12, 14, 16
  };

    const std::vector<unsigned int> timingMValues = {
        1, 2, 5, 10, 20, 30, 50
    };

    const long seed = 20171215;
    const unsigned int numRuns = 11;
/*
    runTimingForMethod<FriedelKellerLibseq>(
        "Friedel-Keller/libseq",
        "friedel-keller-libseq",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);

    runTimingForMethod<BurleyQuantLib>(
        "Burley/QuantLib",
        "burley-quantlib",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);

    runTimingForMethod<PbrtOwenScrambler>(
        "PBRT-OwenScrambler",
        "pbrt-owen",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);

    runTimingForMethod<PbrtFastOwenScrambler>(
        "PBRT-FastOwenScrambler",
        "pbrt-fast-owen",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);

*/
    runTimingForMethod<HelmerStochasticSobol, true>(
        "Helmer-Stochastic-Sobol",
        "helmer-stochastic-sobol",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);
    runTimingForMethod<NewCachedNusPointSet, true>(
        "NewCachedNusPointSet",
        "new-cached-nus",
        models,
        sValues,
        kValues,
        timingMValues,
        seed,
        numRuns,
        outputDir);

    // runTimingForMethod<AdaptiveRegularTiles>(
    // "Adaptive-Regular-Tiles",
    // "adaptive-regular-tiles",
    // models,
    // sValues,
    // kValues,
    // timingMValues,
    // seed,
    // numRuns,
    // outputDir);

    return 0;
}
