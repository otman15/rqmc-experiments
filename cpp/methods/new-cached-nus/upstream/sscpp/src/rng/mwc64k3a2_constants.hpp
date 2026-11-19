/**
 * @file mwc64k3a2_constants.hpp
 * @brief Fixed jump constants of `MWC64k3a2` (internal header).
 *
 * GENERATED FILE, DO NOT EDIT. Produced by tools/mwc64k3a2_constants.py.
 * Parameters: A2 = 184698970548483715, A3 = 6028691832887,
 * m = A3*b^3 + A2*b^2 - 1, b = 2^64, Montgomery radix R = b^4 = 2^256.
 * Each number is stored as 4 little-endian 64-bit words.
 */
#pragma once

#include <array>
#include <cstdint>

namespace sscpp::rng::detail::mwc64k3a2 {

/// Stream jump multiplier in Montgomery form: (b^-1)^(2^169) * R mod m.
inline constexpr std::array<std::uint64_t, 4> kStreamJump = {0xeb1ddd0516d5e95bULL, 0x1b87e1be017d4fa2ULL, 0x8e4285b018082872ULL, 0x0000042274fc84d3ULL};

/// Substream jump multiplier in Montgomery form: (b^-1)^(2^118) * R mod m.
inline constexpr std::array<std::uint64_t, 4> kSubstreamJump = {0xe86e75b5d45672c8ULL, 0xc48c9de952cfed7bULL, 0xa80e44ddb46ba79aULL, 0x000002b23e0a874eULL};

/// R^2 mod m, used to convert a seed to its LCG state.
inline constexpr std::array<std::uint64_t, 4> kR2 = {0x87f44e0a1524463dULL, 0xaf555b19e8c7ac3aULL, 0xd81ad77ec181dc5cULL, 0x00000511a87b638dULL};

/// LCG state y (normal form) of the default package seed {1, 3, 4, 5}.
inline constexpr std::array<std::uint64_t, 4> kDefaultSeedY = {0x0000000000000001ULL, 0x0000000000000003ULL, 0xfd6fd143ce13c981ULL, 0x0000000000000004ULL};

}  // namespace sscpp::rng::detail::mwc64k3a2
