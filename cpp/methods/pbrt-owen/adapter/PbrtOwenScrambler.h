#ifndef PBRT_OWEN_SCRAMBLER_H
#define PBRT_OWEN_SCRAMBLER_H

#include <cstdint>
#include <vector>

class PbrtOwenScrambler
{
public:
    PbrtOwenScrambler(unsigned int s, unsigned int k, long seed);

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    static std::uint32_t scramble(
        std::uint32_t value,
        std::uint32_t seed);

    void updateCurrentPoint();

    unsigned int dimension_;
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;
    long seed_;

    std::vector<double> currentPoint_;
};

#endif
