#ifndef HELMER_STOCHASTIC_SOBOL_H
#define HELMER_STOCHASTIC_SOBOL_H

#include <cstdint>
#include <vector>

class HelmerStochasticSobol
{
public:
    HelmerStochasticSobol(
        unsigned int s,
        unsigned int k,
        long seed);

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    unsigned int dimension_;
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::vector<double> points_;
};

#endif
