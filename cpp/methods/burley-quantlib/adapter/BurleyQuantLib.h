#ifndef BURLEY_QUANTLIB_H
#define BURLEY_QUANTLIB_H

#include <cstdint>
#include <memory>
#include <vector>

namespace QuantLib {
class Burley2020SobolRsg;
}

class BurleyQuantLib
{
public:
    BurleyQuantLib(unsigned int s, unsigned int k, long seed);
    ~BurleyQuantLib();

    BurleyQuantLib(const BurleyQuantLib&) = delete;
    BurleyQuantLib& operator=(const BurleyQuantLib&) = delete;

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::unique_ptr<QuantLib::Burley2020SobolRsg> sequence_;
    std::vector<double> currentPoint_;
};

#endif
