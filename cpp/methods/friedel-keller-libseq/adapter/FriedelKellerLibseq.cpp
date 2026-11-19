#include "FriedelKellerLibseq.h"

#include <cstdlib>

#include "RandomizedTSSequence.h"

FriedelKellerLibseq::FriedelKellerLibseq(
    unsigned int s,
    unsigned int k,
    long seed)
    : numPoints_(std::uint64_t{1} << k),
      currentIndex_(0)
{
    srand48(seed);

    ring_ = std::make_unique<Ring>(2);

    char generator[] = "Niederreiter";
    matrix_ = std::make_unique<C_matrix>(
        ring_.get(), generator, s);

    sequence_ = std::make_unique<RandomizedTSSequence_base_2>(
        matrix_.get(), s, numPoints_);
}

FriedelKellerLibseq::~FriedelKellerLibseq() = default;

void FriedelKellerLibseq::randomize()
{
    sequence_->random_restart();
    currentIndex_ = 0;
}

void FriedelKellerLibseq::reseed(std::uint64_t seed)
{
    srand48(static_cast<long>(static_cast<std::uint32_t>(seed)));
}

void FriedelKellerLibseq::reset()
{
    sequence_->restart();
    currentIndex_ = 0;
}

double FriedelKellerLibseq::coordinate(unsigned int j) const
{
    return (*sequence_)[j];
}

void FriedelKellerLibseq::next()
{
    if (currentIndex_ + 1 < numPoints_)
    {
        ++(*sequence_);
        ++currentIndex_;
    }
}

std::uint64_t FriedelKellerLibseq::numPoints() const
{
    return numPoints_;
}
