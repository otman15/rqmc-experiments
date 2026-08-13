#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BurleyQuantLib.h"
#include "FriedelKellerLibseq.h"
#include "HelmerStochasticSobol.h"
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

template <typename Sequence>
void runStatisticsForMethod(
    const std::string& methodName,
    const std::vector<std::pair<std::string, std::string>>& models,
    const std::vector<unsigned int>& sValues,
    const std::vector<unsigned int>& kValues,
    unsigned int statisticsM,
    long seed,
    const std::filesystem::path& outputDir)
{
    for (const auto& [modelName, model] : models)
    {
        for (unsigned int s : sValues)
        {
            for (unsigned int k : kValues)
            {
                const auto result =
                    runRqmcReplications<Sequence>(
                        model,
                        k,
                        s,
                        statisticsM,
                        seed);

                const double avg =
                    mean(result.estimates);

                const double variance =
                    sampleVariance(result.estimates);

                const double skew =
                    skewness(result.estimates);

                const double kurtosis =
                    excessKurtosis(result.estimates);

                std::cout
                    << "method=" << methodName
                    << ", model=" << modelName
                    << ", s=" << s
                    << ", k=" << k
                    << ", m=" << statisticsM
                    << ", mean="
                    << std::scientific
                    << std::setprecision(16)
                    << avg
                    << ", variance=" << variance
                    << ", skewness=" << skew
                    << ", excess kurtosis=" << kurtosis
                    << ", time="
                    << std::fixed
                    << std::setprecision(6)
                    << result.elapsedTimeSeconds
                    << " seconds\n";

                std::vector<double> sortedEstimates =
                    result.estimates;

                std::sort(
                    sortedEstimates.begin(),
                    sortedEstimates.end());

                const std::filesystem::path outputFile =
                    outputDir /
                    (modelName + "-" +
                     std::to_string(s) + "-" +
                     methodName + "-" +
                     std::to_string(k) + "-" +
                     std::to_string(statisticsM) +
                     ".dat");

                std::ofstream file(outputFile);

                if (!file)
                    throw std::runtime_error(
                        "Cannot open output file: " +
                        outputFile.string());

                file << std::setprecision(17);

                for (double estimate : sortedEstimates)
                    file << estimate << '\n';

                std::cout
                    << "Results written to "
                    << outputFile << '\n';
            }
        }
    }
}

int main()
{
    const std::filesystem::path outputDir =
        "cpp-results/dat";

    std::filesystem::create_directories(outputDir);

  const std::vector<std::pair<std::string, std::string>> models = {
      {"SUMUEU", "sumueu"},
      {"SMOOTHPERB4", "smoothperb4"},
      {"MC2", "mc2"},
      {"POLYNOMIAL", "polynomial"}
  };

  const std::vector<unsigned int> sValues = {
      2, 4, 8, 16, 32
  };

  const std::vector<unsigned int> kValues = {
      8, 10, 12, 14, 16
  };

    const unsigned int statisticsM = 1000;
    const long seed = 20171215;
/*
    runStatisticsForMethod<FriedelKellerLibseq>(
        "Friedel-Keller-libseq",
        models,
        sValues,
        kValues,
        statisticsM,
        seed,
        outputDir);

    runStatisticsForMethod<BurleyQuantLib>(
        "Burley-QuantLib",
        models,
        sValues,
        kValues,
        statisticsM,
        seed,
        outputDir);

    runStatisticsForMethod<PbrtOwenScrambler>(
        "PBRT-OwenScrambler",
        models,
        sValues,
        kValues,
        statisticsM,
        seed,
        outputDir);

    runStatisticsForMethod<PbrtFastOwenScrambler>(
        "PBRT-FastOwenScrambler",
        models,
        sValues,
        kValues,
        statisticsM,
        seed,
        outputDir);

    runStatisticsForMethod<HelmerStochasticSobol>(
        "Helmer-Stochastic-Sobol",
        models,
        sValues,
        kValues,
        statisticsM,
        seed,
        outputDir);
*/        
    runStatisticsForMethod<AdaptiveRegularTiles>(
    "Adaptive-Regular-Tiles",
    models,
    sValues,
    kValues,
    statisticsM,
    seed,
    outputDir);

    return 0;
}
