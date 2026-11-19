#ifndef HELMER_NATIVE_SOBOL_H
#define HELMER_NATIVE_SOBOL_H

#include <cstdint>
#include <vector>

// Timing wrapper around the unchanged upstream GetStochasticSobolSamples().
// Construction allocates the output buffer only. Upstream controls seeding;
// there is deliberately no seed or reseed API.
class HelmerNativeSobol
{
public:
    HelmerNativeSobol(unsigned int s, unsigned int k);

    // Calls upstream with shuffle=false, candidates=1, and owen=false.
    // Includes its RNG initialization and temporary allocations on every call.
    void randomize();
    void reset();
    double coordinate(unsigned int j) const;
    void next();
    std::uint64_t numPoints() const;

private:
    unsigned int dimension_;
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_ = 0;
    std::vector<double> points_;
};

#endif
