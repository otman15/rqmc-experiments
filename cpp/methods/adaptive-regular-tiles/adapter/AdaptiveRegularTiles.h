#ifndef ADAPTIVE_REGULAR_TILES_H
#define ADAPTIVE_REGULAR_TILES_H

#include <cstdint>
#include <vector>

class AdaptiveRegularTiles
{
public:
    AdaptiveRegularTiles(
        unsigned int s,
        unsigned int k,
        long seed);

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    static std::uint32_t scramble(
        std::uint32_t value,
        const std::uint32_t* dna);

    void updateCurrentPoint();

    unsigned int dimension_;
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::vector<std::uint32_t> dna_;
    std::vector<double> currentPoint_;
};

#endif
