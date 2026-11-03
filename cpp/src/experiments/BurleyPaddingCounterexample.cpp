// Standalone, deterministic counterexample to preservation of a global
// (3,5,8)-net under four-dimensional shuffled Sobol padding.
// Compile from the cpp directory:
//   g++ -std=c++17 -O2 src/experiments/BurleyPaddingCounterexample.cpp -o /tmp/burley_padding_counterexample
// Run from the cpp directory so the data path below resolves correctly.
//   /tmp/burley_padding_counterexample             # original fixed seeds
//   /tmp/burley_padding_counterexample 12345       # s0=0, s1=12345
//   /tmp/burley_padding_counterexample 7 12345     # s0=7, s1=12345

// The bundled data/sobol/new-joe-kuo-6.21201 direction-number file is from
// Stephen Joe and Frances Y. Kuo, https://web.maths.unsw.edu.au/~fkuo/sobol/ .
// Their redistribution notice for sobol.cc and its accompanying direction
// numbers follows (the data file itself is kept unchanged):
//
// Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the names of the copyright holders nor the names of the
//       University of New South Wales and the University of Waikato
//       and its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

constexpr const char* JOE_KUO_PATH = "data/sobol/new-joe-kuo-6.21201";
constexpr int DIMENSIONS = 8;
constexpr int BITS = 32;
constexpr int N = 32;  // 2^m, with m = 5.
constexpr std::uint32_t DEFAULT_S0 = 0u;
constexpr std::uint32_t DEFAULT_S1 = 1228959102u;

using Counts = std::array<int, 4>;  // Cells in the order 00, 01, 10, 11.
using Point = std::array<std::uint32_t, DIMENSIONS>;
using PointSet = std::array<Point, N>;
using Directions = std::array<std::array<std::uint32_t, BITS>, DIMENSIONS>;

// An always-enabled assertion: failure must also stop a release build.
void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cout << std::flush;
        std::cerr << "ASSERTION FAILED: " << message << '\n';
        std::abort();
    }
    assert(condition);
}

// Accept decimal 32-bit seeds and reject signs, trailing text, and overflow.
std::uint32_t parseSeed(const char* text, const char* name) {
    check(*text != '\0', std::string(name) + " is empty");
    std::uint64_t value = 0;
    for (const char* digit = text; *digit != '\0'; ++digit) {
        check(*digit >= '0' && *digit <= '9',
              std::string(name) + " must be a decimal uint32_t");
        value = value * 10 + std::uint64_t(*digit - '0');
        check(value <= 0xffffffffu, std::string(name) + " exceeds uint32_t");
    }
    return static_cast<std::uint32_t>(value);
}

struct JoeKuoRow {
    int degree = 0;
    std::uint32_t polynomial = 0;
    std::vector<std::uint32_t> initial;
};

// The official file has one header, then rows d s a m_1 ... m_s.
// Mathematical dimension 1 (C++ dimension 0) is the special radical inverse.
std::array<JoeKuoRow, DIMENSIONS> loadJoeKuo() {
    std::ifstream input(JOE_KUO_PATH);
    check(input.good(), std::string("cannot open Joe-Kuo file: ") + JOE_KUO_PATH);

    std::string line;
    check(static_cast<bool>(std::getline(input, line)), "missing Joe-Kuo header");
    std::array<JoeKuoRow, DIMENSIONS> rows{};
    for (int dim = 1; dim < DIMENSIONS; ++dim) {
        check(static_cast<bool>(std::getline(input, line)),
              "Joe-Kuo file ends before mathematical dimension " + std::to_string(dim + 1));
        std::istringstream fields(line);
        int fileDimension = 0;
        auto& row = rows[dim];
        check(static_cast<bool>(fields >> fileDimension >> row.degree >> row.polynomial),
              "malformed Joe-Kuo row for dimension " + std::to_string(dim + 1));
        check(fileDimension == dim + 1, "unexpected Joe-Kuo dimension number");
        check(row.degree >= 1 && row.degree < BITS, "invalid polynomial degree");
        check(row.polynomial < (std::uint32_t{1} << row.degree),
              "invalid polynomial coefficients");
        for (int column = 0; column < row.degree; ++column) {
            std::uint32_t value = 0;
            check(static_cast<bool>(fields >> value), "missing initial direction number");
            check((value & 1u) != 0u &&
                      value < (std::uint64_t{1} << (column + 1)),
                  "invalid initial direction number");
            row.initial.push_back(value);
        }
        fields >> std::ws;
        check(fields.eof(), "extra fields in Joe-Kuo row");
    }
    return rows;
}

// Column b of C_dim is stored as a uint32_t direction number; its MSB is
// the first output digit. This is the standard Joe-Kuo recurrence, to 32 bits.
Directions buildDirectionNumbers(const std::array<JoeKuoRow, DIMENSIONS>& rows) {
    Directions directions{};
    for (int bit = 0; bit < BITS; ++bit)
        directions[0][bit] = std::uint32_t{1} << (31 - bit);

    for (int dim = 1; dim < DIMENSIONS; ++dim) {
        const auto& row = rows[dim];
        auto& v = directions[dim];
        for (int bit = 0; bit < row.degree; ++bit)
            v[bit] = row.initial[bit] << (31 - bit);
        for (int bit = row.degree; bit < BITS; ++bit) {
            v[bit] = v[bit - row.degree] ^ (v[bit - row.degree] >> row.degree);
            for (int k = 1; k < row.degree; ++k) {
                if ((row.polynomial >> (row.degree - 1 - k)) & 1u)
                    v[bit] ^= v[bit - k];
            }
        }
    }
    return directions;
}

// X_dim(index) = C_dim * binary(index) over F_2. In particular, the index
// bits are NOT Gray coded, and all 32 bits of a transformed index are used.
std::uint32_t sobolCoordinate(const Directions& directions, int dim,
                              std::uint32_t index) {
    std::uint32_t coordinate = 0;
    for (int bit = 0; bit < BITS; ++bit) {
        if ((index >> bit) & 1u)
            coordinate ^= directions[dim][bit];
    }
    return coordinate;
}

void verifyFirstRows(const Directions& directions) {
    constexpr std::array<const char*, DIMENSIONS> expected = {
        "10000", "11111", "11011", "11001",
        "10011", "10001", "11110", "10101"
    };
    std::cout << "Joe-Kuo first-row convention check (a0...a4):\n";
    for (int dim = 0; dim < DIMENSIONS; ++dim) {
        std::string form;
        for (int bit = 0; bit < 5; ++bit)
            form += char('0' + ((directions[dim][bit] >> 31) & 1u));
        std::cout << "  dim " << dim << " : " << form << '\n';
        check(form == expected[dim],
              "Joe-Kuo data or Gray-code/index convention is wrong at dim " +
                  std::to_string(dim));
    }
}

std::string formatCounts(const Counts& counts) {
    std::ostringstream out;
    out << '[' << counts[0] << ',' << counts[1] << ','
        << counts[2] << ',' << counts[3] << ']';
    return out.str();
}

Counts countQuadrants(const PointSet& points, int firstDim, int secondDim) {
    Counts counts{};
    for (const auto& point : points) {
        const unsigned first = (point[firstDim] >> 31) & 1u;
        const unsigned second = (point[secondDim] >> 31) & 1u;
        ++counts[(first << 1) | second];
    }
    return counts;
}

struct NetChecks {
    int quarterPasses = 0;
    int quadrantPasses = 0;

    bool isNet() const { return quarterPasses == 8 && quadrantPasses == 28; }
};

// Check every volume-1/4 elementary interval of a transformed 8D set.
// A supplied seed pair may pass or fail; report the result without presuming it.
NetChecks inspectNet(const PointSet& points, const char* name) {
    constexpr Counts uniform{8, 8, 8, 8};
    NetChecks result;
    std::cout << "\n=== Full (3,5,8)-net check: " << name << " ===\n";
    for (int dim = 0; dim < DIMENSIONS; ++dim) {
        Counts counts{};
        for (const auto& point : points)
            ++counts[(point[dim] >> 30) & 3u];
        if (counts == uniform) {
            ++result.quarterPasses;
        } else {
            std::cout << "  dim " << dim << " quarters: "
                      << formatCounts(counts) << " FAIL\n";
        }
    }
    for (int first = 0; first < DIMENSIONS; ++first) {
        for (int second = first + 1; second < DIMENSIONS; ++second) {
            const Counts counts = countQuadrants(points, first, second);
            if (counts == uniform) {
                ++result.quadrantPasses;
            } else {
                std::cout << "  dims " << first << ',' << second
                          << " quadrants: " << formatCounts(counts) << " FAIL\n";
            }
        }
    }
    std::cout << "one-coordinate quarter tests: " << result.quarterPasses << "/8 PASS\n"
              << "two-coordinate quadrant tests: " << result.quadrantPasses << "/28 PASS\n"
              << "==> " << (result.isNet() ? "satisfies" : "fails")
              << " the (3,5,8)-net condition\n";
    return result;
}

// For m-t = 2, these eight quarter and 28 quadrant checks enumerate every
// base-two elementary interval of volume 2^-2 in eight dimensions.
void verifyReferenceNet(const PointSet& reference) {
    constexpr Counts uniform{8, 8, 8, 8};
    int quarterPasses = 0;
    int quadrantPasses = 0;
    std::cout << "\n=== Reference 8D Joe-Kuo net ===\n";
    for (int dim = 0; dim < DIMENSIONS; ++dim) {
        Counts counts{};
        for (const auto& point : reference)
            ++counts[(point[dim] >> 30) & 3u];
        std::cout << "  dim " << dim << " quarters: " << formatCounts(counts) << '\n';
        check(counts == uniform, "reference quarter test failed at dim " +
                                     std::to_string(dim));
        ++quarterPasses;
    }
    for (int first = 0; first < DIMENSIONS; ++first) {
        for (int second = first + 1; second < DIMENSIONS; ++second) {
            const Counts counts = countQuadrants(reference, first, second);
            std::cout << "  dims " << first << ',' << second
                      << " quadrants: " << formatCounts(counts) << '\n';
            check(counts == uniform,
                  "reference quadrant test failed at dims " +
                      std::to_string(first) + "," + std::to_string(second));
            ++quadrantPasses;
        }
    }
    check(quarterPasses == 8 && quadrantPasses == 28,
          "wrong number of reference net tests");
    std::cout << "one-coordinate quarter tests: " << quarterPasses << "/8 PASS\n"
              << "two-coordinate quadrant tests: " << quadrantPasses << "/28 PASS\n"
              << "==> verified (3,5,8)-net\n";
}

std::uint32_t hash_combine(std::uint32_t seed, std::uint32_t value) {
    return seed ^ (value + (seed << 6) + (seed >> 2));
}

std::uint32_t reverse_bits(std::uint32_t x) {
    x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
    x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
    x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
    x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
    return (x >> 16) | (x << 16);
}

std::uint32_t laine_karras_permutation(std::uint32_t x, std::uint32_t seed) {
    x += seed;
    x ^= x * 0x6c50b47cu;
    x ^= x * 0xb82f1e52u;
    x ^= x * 0xc7afe638u;
    x ^= x * 0x8d22f6e6u;
    return x;
}

// P_s on an index is a nested shuffle; on a coordinate it is a nested scramble.
std::uint32_t nested_uniform_scramble_base2(std::uint32_t x, std::uint32_t seed) {
    x = reverse_bits(x);
    x = laine_karras_permutation(x, seed);
    x = reverse_bits(x);
    return x;
}

std::uint32_t verifyAlignedBlock(const std::array<std::uint32_t, N>& indices,
                                 const char* label) {
    const std::uint32_t block = indices[0] >> 5;
    std::array<bool, N> seen{};
    for (const auto index : indices) {
        check((index >> 5) == block,
              std::string(label) + " does not stay in one aligned block");
        const std::uint32_t low = index & 31u;
        check(!seen[low], std::string(label) + " repeats a low-five-bit index");
        seen[low] = true;
    }
    for (bool wasSeen : seen)
        check(wasSeen, std::string(label) + " omits a low-five-bit index");
    std::cout << label << " aligned block q = " << block
              << ", low five bits are a permutation of 0...31\n";
    return block;
}

std::string firstBitPair(const Point& point, int firstDim, int secondDim) {
    std::string result = "00";
    result[0] += char((point[firstDim] >> 31) & 1u);
    result[1] += char((point[secondDim] >> 31) & 1u);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc > 3) {
        std::cerr << "Usage: " << argv[0] << " [s1 | s0 s1]\n";
        return 2;
    }
    const std::uint32_t s0 = argc == 3 ? parseSeed(argv[1], "s0") : DEFAULT_S0;
    const std::uint32_t s1 = argc == 3 ? parseSeed(argv[2], "s1") :
                             argc == 2 ? parseSeed(argv[1], "s1") : DEFAULT_S1;
    const bool fixedCounterexample = s0 == DEFAULT_S0 && s1 == DEFAULT_S1;

    const auto rows = loadJoeKuo();
    const Directions directions = buildDirectionNumbers(rows);
    std::cout << "Joe-Kuo file: " << JOE_KUO_PATH << '\n';
    verifyFirstRows(directions);

    PointSet reference{};
    for (std::uint32_t i = 0; i < N; ++i) {
        for (int dim = 0; dim < DIMENSIONS; ++dim)
            reference[i][dim] = sobolCoordinate(directions, dim, i);
    }
    verifyReferenceNet(reference);  // Stop before any Burley operation on failure.

    std::cout << "\nGroup seeds: s0 = " << s0 << ", s1 = " << s1 << '\n';
    std::array<std::uint32_t, N> p0{};
    std::array<std::uint32_t, N> p1{};
    for (std::uint32_t i = 0; i < N; ++i) {
        p0[i] = nested_uniform_scramble_base2(i, s0);
        p1[i] = nested_uniform_scramble_base2(i, s1);
    }
    std::cout << "\n=== Hierarchical index-block checks ===\n";
    verifyAlignedBlock(p0, "s0");
    verifyAlignedBlock(p1, "s1");

    const std::uint32_t coordinateSeed0 = hash_combine(s0, 3u);
    const std::uint32_t coordinateSeed1 = hash_combine(s1, 3u);
    std::cout << "\ncoordinate seeds for local dimension 3:\n"
              << "  hash_combine(s0,3) = " << coordinateSeed0 << '\n'
              << "  hash_combine(s1,3) = 0x" << std::hex << coordinateSeed1
              << std::dec << '\n';
    if (fixedCounterexample) {
        check(coordinateSeed0 == 3u, "unexpected coordinate seed for group 0");
        check(coordinateSeed1 == 0x2b2b179cu,
              "unexpected coordinate seed for group 1");
    }

    PointSet controlled{};
    PointSet literalBefore{};
    PointSet literalAfter{};
    for (std::uint32_t i = 0; i < N; ++i) {
        for (int local = 0; local < 4; ++local) {
            const std::uint32_t first = sobolCoordinate(directions, local, p0[i]);
            const std::uint32_t second = sobolCoordinate(directions, local, p1[i]);

            // Controlled 8D set: retain C0...C7, but use P_s0 and P_s1
            // for the respective four-dimensional blocks.
            controlled[i][local] = first;
            controlled[i][local + 4] =
                sobolCoordinate(directions, local + 4, p1[i]);

            // Literal padding: reuse C0...C3 in the second block.
            literalBefore[i][local] = first;
            literalBefore[i][local + 4] = second;

            // Apply the specified coordinate scramble to every coordinate.
            literalAfter[i][local] = nested_uniform_scramble_base2(
                first, hash_combine(s0, std::uint32_t(local)));
            literalAfter[i][local + 4] = nested_uniform_scramble_base2(
                second, hash_combine(s1, std::uint32_t(local)));
        }
    }

    const Counts controlledCounts = countQuadrants(controlled, 3, 7);
    const Counts beforeCounts = countQuadrants(literalBefore, 3, 7);
    const Counts afterCounts = countQuadrants(literalAfter, 3, 7);
    constexpr Counts uniform{8, 8, 8, 8};
    if (fixedCounterexample) {
        check(controlledCounts == Counts{9, 7, 7, 9},
              "controlled counts differ from the fixed-seed reference: " +
                  formatCounts(controlledCounts));
        check(beforeCounts == Counts{7, 9, 9, 7},
              "literal-padding pre-scramble counts differ: " +
                  formatCounts(beforeCounts));
        check(afterCounts == Counts{9, 7, 7, 9},
              "literal-padding post-scramble counts differ: " +
                  formatCounts(afterCounts));
        check(controlledCounts != uniform && beforeCounts != uniform &&
                  afterCounts != uniform,
              "a fixed-seed failing cross-block projection is actually uniform");
    }

    std::cout << "\n=== Computed cross-block quadrant counts [00,01,10,11] ===\n"
              << "Expected for a (3,5,8)-net: " << formatCounts(uniform) << '\n'
              << "Controlled different index maps (dims 3,7): "
              << formatCounts(controlledCounts) << '\n'
              << "Literal padding before scramble (dims 3,7): "
              << formatCounts(beforeCounts) << '\n'
              << "Literal padding after scramble (dims 3,7): "
              << formatCounts(afterCounts) << '\n';

    const NetChecks controlledNet = inspectNet(controlled, "controlled index maps");
    const NetChecks beforeNet = inspectNet(literalBefore, "literal padding before scramble");
    const NetChecks afterNet = inspectNet(literalAfter, "literal padding after scramble");
    if (fixedCounterexample)
        check(!afterNet.isNet(), "fixed-seed padded set unexpectedly passed all net tests");

    std::cout << "\n=== Reproducibility table ===\n"
              << " i         P_s0(i)      P_s1(i)  low0 low1 controlled before after\n";
    for (int i = 0; i < N; ++i) {
        std::cout << std::setw(2) << i << ' '
                  << std::setw(15) << p0[i] << ' '
                  << std::setw(12) << p1[i] << ' '
                  << std::setw(5) << (p0[i] & 31u) << ' '
                  << std::setw(4) << (p1[i] & 31u) << ' '
                  << std::setw(10) << firstBitPair(controlled[i], 3, 7) << ' '
                  << std::setw(6) << firstBitPair(literalBefore[i], 3, 7) << ' '
                  << firstBitPair(literalAfter[i], 3, 7) << '\n';
    }

    std::cout << "\n=== Final summary ===\n"
              << "Reference Joe-Kuo 8D: 8/8 quarter and 28/28 quadrant tests PASS; "
                 "(3,5,8)-net.\n"
              << "Group seeds: s0 = " << s0 << ", s1 = " << s1 << '\n'
              << "Controlled different-index-map counts: "
              << formatCounts(controlledCounts) << "; full net: "
              << (controlledNet.isNet() ? "PASS" : "FAIL") << ".\n"
              << "Literal padding before coordinate scramble: "
              << formatCounts(beforeCounts) << "; full net: "
              << (beforeNet.isNet() ? "PASS" : "FAIL") << ".\n"
              << "Literal padding after coordinate scramble: "
              << formatCounts(afterCounts) << "; full net: "
              << (afterNet.isNet() ? "PASS" : "FAIL") << ".\n";
    if (fixedCounterexample)
        std::cout << "Fixed-seed counterexample reproduced successfully.\n";
    else if (afterNet.isNet())
        std::cout << "For these seeds, the padded set passes all (3,5,8)-net tests.\n";
    else
        std::cout << "For these seeds, Burley-style padding fails the global "
                     "(3,5,8)-net condition.\n";
}
