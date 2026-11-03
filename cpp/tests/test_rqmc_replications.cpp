#include "RqmcReplications.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "HelmerStochasticSobol.h"
#include "NewCachedNusPointSet.h"

namespace
{
// Release builds must execute all checks.
void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

double cpuNow()
{
    timespec time;
    require(clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time) == 0, "process clock failed");
    return static_cast<double>(time.tv_sec) + static_cast<double>(time.tv_nsec) * 1e-9;
}

void burnCpu()
{
    // Non-elidable work makes missing phase time (e.g. just the final
    // randomization) detectable without a machine-dependent speed threshold.
    volatile std::uint64_t value = 1;
    for (unsigned int i = 0; i < 20000; ++i)
        value = value * 6364136223846793005ULL + 1;
}

struct Counts
{
    unsigned int factories = 0, constructions = 0, destructions = 0;
    unsigned int randomizations = 0, resets = 0, advances = 0, coordinates = 0;
    double preparationInside = 0.0, scramblingInside = 0.0, evaluationInside = 0.0;
} counts;

double probeCoordinate(unsigned int rep, std::uint64_t i, unsigned int j)
{
    return static_cast<double>((11 + rep * 17 + i * 3 + j) % 64) / 64;
}

class ProbeSequence
{
public:
    ProbeSequence() : points_(12)
    {
        const auto start = cpuNow();
        ++counts.constructions;
        burnCpu();
        counts.preparationInside += cpuNow() - start;
    }
    ~ProbeSequence() { ++counts.destructions; }
    ProbeSequence(const ProbeSequence&) = delete;
    ProbeSequence(ProbeSequence&&) = delete;
    std::uint64_t numPoints() const { return 4; }
    void randomize()
    {
        const auto start = cpuNow();
        ++counts.randomizations;
        ++rep_;
        burnCpu();
        for (std::uint64_t i = 0; i < 4; ++i)
            for (unsigned int j = 0; j < 3; ++j)
                points_[i * 3 + j] = probeCoordinate(rep_, i, j);
        current_ = 3; // Deliberately does not rewind: the runner must reset.
        counts.scramblingInside += cpuNow() - start;
    }
    void reset() { ++counts.resets; current_ = 0; }
    void next() { ++counts.advances; ++current_; }
    double coordinate(unsigned int j) const
    {
        const auto start = cpuNow();
        require(current_ < 4 && j < 3, "runner accessed outside the current point");
        ++counts.coordinates;
        burnCpu();
        const auto value = points_[current_ * 3 + j];
        counts.evaluationInside += cpuNow() - start;
        return value;
    }
private:
    std::vector<double> points_;
    unsigned int current_ = 0, rep_ = 0;
};

void checkTimings(const RqmcReplicationResult& result)
{
    const double phases = result.preparationTimeSeconds + result.scramblingTimeSeconds +
                          result.evaluationTimeSeconds;
    require(std::isfinite(result.elapsedTimeSeconds) &&
            result.elapsedTimeSeconds >= 0 && result.preparationTimeSeconds >= 0 &&
            result.scramblingTimeSeconds >= 0 && result.evaluationTimeSeconds >= 0,
            "invalid CPU timings");
    require(result.elapsedTimeSeconds + 1e-8 >= phases, "complete experiment omitted a phase");
}

void testPreparationAndTiming(unsigned int m)
{
    counts = {};
    const auto result = runPreparedRqmcReplications("mc2", 3, m, []
    {
        ++counts.factories;
        return ProbeSequence(); // Nonmovable, like the cached NUS adapter.
    });
    require(counts.factories == 1 && counts.constructions == 1 && counts.destructions == 1,
            "prepare/construct/destroy must happen exactly once per run");
    require(counts.randomizations == m, "must randomize exactly once per measured replication");
    require(counts.resets == m, "must rewind after every randomization");
    require(counts.advances == m * 3 && counts.coordinates == m * 12,
            "must evaluate every point and dimension exactly once");
    require(result.estimates.size() == m, "wrong replication estimate count");
    std::vector<double> x(3);
    for (unsigned int rep = 0; rep < m; ++rep)
    {
        double sum = 0;
        for (std::uint64_t i = 0; i < 4; ++i)
        {
            for (unsigned int j = 0; j < 3; ++j)
                x[j] = probeCoordinate(rep + 1, i, j);
            sum += evaluateRqmcModel("mc2", x);
        }
        require(result.estimates[rep] == sum / 4, "incorrect estimate or iteration order");
    }
    checkTimings(result);
    // The inner timers exclude runner overhead and must fit entirely inside
    // their reported phase totals, for every call, not just the last call.
    require(result.preparationTimeSeconds + 1e-8 >= counts.preparationInside, "preparation omitted work");
    require(result.scramblingTimeSeconds + 1e-8 >= counts.scramblingInside, "scrambling omitted randomizations");
    require(result.evaluationTimeSeconds + 1e-8 >= counts.evaluationInside, "evaluation omitted replications");
    if (m == 0)
        require(result.scramblingTimeSeconds == 0 && result.evaluationTimeSeconds == 0,
                "zero replications must not scramble/evaluate");
}

template <typename Prepare>
std::vector<double> manualEstimates(const std::string& model, unsigned int s, unsigned int m, Prepare prepare)
{
    auto sequence = prepare();
    std::vector<double> result, x(s);
    for (unsigned int r = 0; r < m; ++r)
    {
        sequence.randomize();
        sequence.reset();
        double sum = 0;
        for (std::uint64_t i = 0; i < sequence.numPoints(); ++i)
        {
            for (unsigned int j = 0; j < s; ++j)
                x[j] = sequence.coordinate(j);
            sum += evaluateRqmcModel(model, x);
            if (i + 1 < sequence.numPoints())
                sequence.next();
        }
        result.push_back(sum / static_cast<double>(sequence.numPoints()));
    }
    return result;
}

void testActualAdapters()
{
    for (const auto k : {0u, 1u, 4u})
        for (const auto m : {0u, 1u, 5u})
            for (const long seed : {0L, -7L})
                for (const std::string model : {"mc2", "gaussian", "polynomial"})
                {
                    const auto prepareHck = [=]
                    {
                        return HelmerStochasticSobol(3, k, seed,
                            HelmerStochasticSobol::Initialization::PrepareOnly);
                    };
                    const auto prepareNus = [=] { return NewCachedNusPointSet(3, k, seed); };
                    const auto hck = runPreparedRqmcReplications(model, 3, m, prepareHck);
                    const auto nus = runPreparedRqmcReplications(model, 3, m, prepareNus);
                    require(hck.estimates == manualEstimates(model, 3, m, prepareHck),
                            "HCK runner differs from manual continuous-RNG lifecycle");
                    require(nus.estimates == manualEstimates(model, 3, m, prepareNus),
                            "cached NUS runner differs from manual continuous-RNG lifecycle");
                    checkTimings(hck);
                    checkTimings(nus);
                }
}

std::vector<long> legacySeeds;
class LegacySequence
{
public:
    LegacySequence(unsigned int, unsigned int, long seed) : seed_(seed)
    {
        legacySeeds.push_back(seed);
    }
    // Deliberately has no randomize/reset API: old adapters must still compile.
    std::uint64_t numPoints() const { return 1; }
    void next() {}
    double coordinate(unsigned int j) const { return (seed_ + j) / 64.0; }
private:
    long seed_;
};

void testLegacyRunner()
{
    legacySeeds.clear();
    const auto result = runRqmcReplications<LegacySequence>("mc2", 0, 2, 4, 11);
    require(legacySeeds == std::vector<long>({11, 12, 13, 14}), "legacy seed+r behavior changed");
    for (unsigned int r = 0; r < 4; ++r)
        require(result.estimates[r] == evaluateRqmcModel("mc2", {(11 + r) / 64.0, (12 + r) / 64.0}),
                "legacy estimate behavior changed");
    require(result.preparationTimeSeconds == 0 && result.scramblingTimeSeconds == 0 &&
            result.evaluationTimeSeconds == 0, "legacy timings should not claim phase separation");
}
}

int main()
{
    try
    {
        for (const auto m : {0u, 1u, 4u})
            testPreparationAndTiming(m);
        testActualAdapters();
        testLegacyRunner();
        std::cout << "RQMC replications: all lifecycle and timing checks passed\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "RQMC replication test failed: " << error.what() << '\n';
        return 1;
    }
}
