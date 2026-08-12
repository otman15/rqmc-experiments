#ifndef FRIEDEL_KELLER_LIBSEQ_H
#define FRIEDEL_KELLER_LIBSEQ_H

#include <cstdint>
#include <memory>

class Ring;
class C_matrix;
class RandomizedTSSequence_base_2;

class FriedelKellerLibseq
{
public:
    FriedelKellerLibseq(unsigned int s, unsigned int k, long seed);
    ~FriedelKellerLibseq();

    FriedelKellerLibseq(const FriedelKellerLibseq&) = delete;
    FriedelKellerLibseq& operator=(const FriedelKellerLibseq&) = delete;

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::unique_ptr<Ring> ring_;
    std::unique_ptr<C_matrix> matrix_;
    std::unique_ptr<RandomizedTSSequence_base_2> sequence_;
};

#endif
