#ifndef FRIEDEL_KELLER_LIBSEQ_H
#define FRIEDEL_KELLER_LIBSEQ_H

#include <cstdint>
#include <memory>
#include <vector>

class Ring;
class C_matrix;
class RandomizedTSSequence_base_2;

class FriedelKellerLibseq
{
public:
    FriedelKellerLibseq(unsigned int s, unsigned int k, long seed);
    // Supply 32-bit generating columns in binary-index order, dimension major:
    // binaryColumns[j * k + c]. Construction still includes the initial scramble.
    FriedelKellerLibseq(unsigned int s, unsigned int k, long seed,
                       const std::vector<std::uint32_t>& binaryColumns);
    ~FriedelKellerLibseq();

    FriedelKellerLibseq(const FriedelKellerLibseq&) = delete;
    FriedelKellerLibseq& operator=(const FriedelKellerLibseq&) = delete;

    // Calls upstream random_restart(), reusing its arrays and rewinding.
    // Upstream construction already generates an initial scrambling.
    void randomize();
    // Affects future randomizations only. Upstream uses the process-global
    // drand48 RNG and the low 32 seed bits; instances share that RNG state.
    void reseed(std::uint64_t seed);
    void reset(); // Rewinds the current scrambling without consuming randomness.

    double coordinate(unsigned int j) const;
    void next();

    std::uint64_t numPoints() const;

private:
    FriedelKellerLibseq(unsigned int s, unsigned int k, long seed,
                       const std::vector<std::uint32_t>* binaryColumns);

    std::uint64_t numPoints_;
    std::uint64_t currentIndex_;

    std::unique_ptr<Ring> ring_;
    std::unique_ptr<C_matrix> matrix_;
    std::unique_ptr<RandomizedTSSequence_base_2> sequence_;
};

#endif
