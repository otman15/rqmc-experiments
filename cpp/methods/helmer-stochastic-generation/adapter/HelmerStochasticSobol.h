#ifndef HELMER_STOCHASTIC_SOBOL_H
#define HELMER_STOCHASTIC_SOBOL_H

#include <cstdint>
#include <memory>
#include <vector>

class HelmerStochasticSobol
{
public:
    enum class Initialization { Randomize, PrepareOnly };

    // Retains the original seeded initial scrambling.
    HelmerStochasticSobol(
        unsigned int s,
        unsigned int k,
        long seed);

    // PrepareOnly allocates a zero cache and seeds the RNG, consuming no
    // scrambling. Call randomize() before evaluating this cache.
    HelmerStochasticSobol(
        unsigned int s,
        unsigned int k,
        long seed,
        Initialization initialization);
    ~HelmerStochasticSobol();
    HelmerStochasticSobol(const HelmerStochasticSobol& other);
    HelmerStochasticSobol& operator=(const HelmerStochasticSobol& other);

    // Reuses storage, advances the upstream PCG32 stream, and resets iteration.
    void randomize();
    // Affects future randomizations, preserving the cache and iteration position.
    void reseed(std::uint64_t seed);
    void reset();

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    friend struct HelmerStochasticSobolTestAccess;
    struct RandomState;
    unsigned int dimension_;
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::vector<double> points_;
    std::unique_ptr<RandomState> random_;
};

#endif
