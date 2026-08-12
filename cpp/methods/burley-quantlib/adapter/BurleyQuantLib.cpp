#include "BurleyQuantLib.h"

#include <ql/math/randomnumbers/burley2020sobolrsg.hpp>

BurleyQuantLib::BurleyQuantLib(
    unsigned int s,
    unsigned int k,
    long seed)
    : numPoints_(std::uint64_t{1} << k),
      currentIndex_(0)
{
    const unsigned long sobolSeed = 42;
    const unsigned long scrambleSeed =
        static_cast<unsigned long>(seed);

    sequence_ =
        std::make_unique<QuantLib::Burley2020SobolRsg>(
            s,
            sobolSeed,
            QuantLib::SobolRsg::Jaeckel,
            scrambleSeed);

    currentPoint_ = sequence_->nextSequence().value;
}

BurleyQuantLib::~BurleyQuantLib() = default;

double BurleyQuantLib::coordinate(unsigned int j) const
{
    return currentPoint_[j];
}

void BurleyQuantLib::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        currentPoint_ = sequence_->nextSequence().value;
        ++currentIndex_;
    }
}

std::uint64_t BurleyQuantLib::numPoints() const
{
    return numPoints_;
}
