#include <cstdint>

#include "rng.h"

namespace
{
// Only compare_nus_methods links these wrappers. The other executables retain
// their native drand48 calls, and no upstream source is modified.
sampling::RNG rng(0);
}

extern "C" double __wrap_drand48()
{
    // Preserve the libraries' double interface using exactly one PCG32 draw.
    return rng.GetUniformFloat();
}

extern "C" void __wrap_srand48(long seed)
{
    // Match the adapters' existing low-32-bit seeding convention. The comparison
    // restricts seeds to this range for Helmer and NewCachedNusPointSet as well.
    rng = sampling::RNG(static_cast<std::uint32_t>(seed));
}
