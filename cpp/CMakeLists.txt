cmake_minimum_required(VERSION 3.16)

project(rqmc_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

set(LIBSEQ_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/friedel-keller-libseq/upstream/libseq_beta_04.21.01)

set(FRIEDEL_KELLER_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/friedel-keller-libseq/adapter)

# ---------------------------------------------------------------------------
# Common experiment code
# ---------------------------------------------------------------------------

add_library(rqmc_common STATIC
    src/common/Integrands.cpp
    src/common/RqmcReplications.cpp
    src/common/Statistics.cpp
)

target_include_directories(rqmc_common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/common
)

# ---------------------------------------------------------------------------
# Original Friedel-Keller libseq
# ---------------------------------------------------------------------------

add_library(libseq STATIC
    ${LIBSEQ_DIR}/own_types.C
    ${LIBSEQ_DIR}/digit_gen.C
    ${LIBSEQ_DIR}/LatinSupercube.C
    ${LIBSEQ_DIR}/LatinHypercube.C
    ${LIBSEQ_DIR}/RandomizedTSSequence_base_2.C
    ${LIBSEQ_DIR}/RandomizedTSSequence.C
    ${LIBSEQ_DIR}/HaltonSequence.C
    ${LIBSEQ_DIR}/DigitalSequence.C
    ${LIBSEQ_DIR}/DigitalSequence_base_2.C
    ${LIBSEQ_DIR}/DigitalSequence_classic.C
    ${LIBSEQ_DIR}/DigitalSequence_medium_base.C
    ${LIBSEQ_DIR}/DigitalSequence_advanced.C
    ${LIBSEQ_DIR}/Sequence.C
    ${LIBSEQ_DIR}/Polynomial.C
    ${LIBSEQ_DIR}/Permutation.C
    ${LIBSEQ_DIR}/Counter.C
    ${LIBSEQ_DIR}/VectorMatrix.C
    ${LIBSEQ_DIR}/VectorRing.C
    ${LIBSEQ_DIR}/C_matrix.C
    ${LIBSEQ_DIR}/Ring.C
)

target_include_directories(libseq PRIVATE
    ${LIBSEQ_DIR}
    ${FRIEDEL_KELLER_ADAPTER_DIR}/compat
)

target_compile_options(libseq PRIVATE
    -O2
    -funroll-loops
)

# ---------------------------------------------------------------------------
# Friedel-Keller adapter
# ---------------------------------------------------------------------------

add_library(friedel_keller_libseq STATIC
    ${FRIEDEL_KELLER_ADAPTER_DIR}/FriedelKellerLibseq.cpp
)

target_include_directories(friedel_keller_libseq PUBLIC
    ${FRIEDEL_KELLER_ADAPTER_DIR}
)

target_include_directories(friedel_keller_libseq PRIVATE
    ${LIBSEQ_DIR}
    ${FRIEDEL_KELLER_ADAPTER_DIR}/compat
)

target_link_libraries(friedel_keller_libseq PRIVATE
    libseq
)

# ---------------------------------------------------------------------------
# Adapter test
# ---------------------------------------------------------------------------

add_executable(test_friedel_keller
    ${FRIEDEL_KELLER_ADAPTER_DIR}/test_adapter.cpp
)

target_link_libraries(test_friedel_keller PRIVATE
    friedel_keller_libseq
)

# ---------------------------------------------------------------------------
# QuantLib - Burley 2020 implementation
# ---------------------------------------------------------------------------

set(QUANTLIB_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/burley-quantlib/upstream/QuantLib-1.43)

set(QL_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(QL_BUILD_TEST_SUITE OFF CACHE BOOL "" FORCE)
set(QL_BUILD_FUZZ_TEST_SUITE OFF CACHE BOOL "" FORCE)

add_subdirectory(
    ${QUANTLIB_DIR}
    ${CMAKE_BINARY_DIR}/quantlib
)

# ---------------------------------------------------------------------------
# Timing experiments
# ---------------------------------------------------------------------------

add_executable(time_experiments
    src/experiments/TimeExperiments.cpp
)

target_link_libraries(time_experiments PRIVATE
    rqmc_common
    friedel_keller_libseq
)

# ---------------------------------------------------------------------------
# Statistics experiments
# ---------------------------------------------------------------------------

add_executable(statistics_experiments
    src/experiments/StatisticsExperiments.cpp
)

target_link_libraries(statistics_experiments PRIVATE
    rqmc_common
    friedel_keller_libseq
)

# ---------------------------------------------------------------------------
# Statistics test
# ---------------------------------------------------------------------------

add_executable(test_statistics
    tests/test_statistics.cpp
)

target_link_libraries(test_statistics PRIVATE
    rqmc_common
)

add_test(NAME statistics COMMAND test_statistics)
add_test(NAME friedel_keller COMMAND test_friedel_keller)
add_test(NAME burley_quantlib COMMAND test_burley_quantlib)
add_test(NAME pbrt_owen COMMAND test_pbrt_owen)
add_test(NAME pbrt_fast_owen COMMAND test_pbrt_fast_owen)
add_test(NAME helmer_stochastic_sobol COMMAND test_helmer_stochastic_sobol)
add_test(NAME adaptive_regular_tiles COMMAND test_adaptive_regular_tiles)

# ---------------------------------------------------------------------------
# QuantLib Burley adapter
# ---------------------------------------------------------------------------

set(BURLEY_QUANTLIB_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/burley-quantlib/adapter)

add_library(burley_quantlib STATIC
    ${BURLEY_QUANTLIB_ADAPTER_DIR}/BurleyQuantLib.cpp
)

target_include_directories(burley_quantlib PUBLIC
    ${BURLEY_QUANTLIB_ADAPTER_DIR}
)

target_link_libraries(burley_quantlib PRIVATE
    ql_library
)

add_executable(test_burley_quantlib
    ${BURLEY_QUANTLIB_ADAPTER_DIR}/test_adapter.cpp
)

target_link_libraries(test_burley_quantlib PRIVATE
    burley_quantlib
)

target_link_libraries(time_experiments PRIVATE
    burley_quantlib
)

target_link_libraries(statistics_experiments PRIVATE
    burley_quantlib
)

# ---------------------------------------------------------------------------
# PBRT Sobol core
# ---------------------------------------------------------------------------

set(PBRT_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/pbrt-owen/upstream/pbrt-v4)

set(PBRT_OWEN_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/pbrt-owen/adapter)

add_library(pbrt_sobol_core STATIC
    ${PBRT_OWEN_ADAPTER_DIR}/PbrtSobolCore.cpp
    ${PBRT_DIR}/src/pbrt/util/sobolmatrices.cpp
)

target_include_directories(pbrt_sobol_core PUBLIC
    ${PBRT_OWEN_ADAPTER_DIR}
    ${PBRT_DIR}/src
)

# ---------------------------------------------------------------------------
# PBRT full Owen scrambler adapter
# ---------------------------------------------------------------------------

add_library(pbrt_owen_scrambler STATIC
    ${PBRT_OWEN_ADAPTER_DIR}/PbrtOwenScrambler.cpp
)

target_include_directories(pbrt_owen_scrambler PUBLIC
    ${PBRT_OWEN_ADAPTER_DIR}
)

target_link_libraries(pbrt_owen_scrambler PRIVATE
    pbrt_sobol_core
)

add_executable(test_pbrt_owen
    ${PBRT_OWEN_ADAPTER_DIR}/test_owen_adapter.cpp
)

target_link_libraries(test_pbrt_owen PRIVATE
    pbrt_owen_scrambler
)

# ---------------------------------------------------------------------------
# PBRT Fast Owen scrambler adapter
# ---------------------------------------------------------------------------

add_library(pbrt_fast_owen_scrambler STATIC
    ${PBRT_OWEN_ADAPTER_DIR}/PbrtFastOwenScrambler.cpp
)

target_include_directories(pbrt_fast_owen_scrambler PUBLIC
    ${PBRT_OWEN_ADAPTER_DIR}
)

target_link_libraries(pbrt_fast_owen_scrambler PRIVATE
    pbrt_sobol_core
)

add_executable(test_pbrt_fast_owen
    ${PBRT_OWEN_ADAPTER_DIR}/test_fast_owen_adapter.cpp
)

target_link_libraries(test_pbrt_fast_owen PRIVATE
    pbrt_fast_owen_scrambler
)

target_link_libraries(time_experiments PRIVATE
    pbrt_owen_scrambler
    pbrt_fast_owen_scrambler
)

target_link_libraries(statistics_experiments PRIVATE
    pbrt_owen_scrambler
    pbrt_fast_owen_scrambler
)

# ---------------------------------------------------------------------------
# Helmer stochastic Sobol adapter
# ---------------------------------------------------------------------------

set(HELMER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/helmer-stochastic-generation/upstream/stochastic-generation)

set(HELMER_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/helmer-stochastic-generation/adapter)

add_library(helmer_stochastic_sobol STATIC
    ${HELMER_ADAPTER_DIR}/HelmerStochasticSobol.cpp
)

target_include_directories(helmer_stochastic_sobol PUBLIC
    ${HELMER_ADAPTER_DIR}
)

target_include_directories(helmer_stochastic_sobol PRIVATE
    ${HELMER_DIR}/sampling
)

add_executable(test_helmer_stochastic_sobol
    ${HELMER_ADAPTER_DIR}/test_adapter.cpp
)

target_link_libraries(test_helmer_stochastic_sobol PRIVATE
    helmer_stochastic_sobol
)

target_include_directories(test_helmer_stochastic_sobol PRIVATE
    ${HELMER_DIR}/sampling
)

target_link_libraries(time_experiments PRIVATE
    helmer_stochastic_sobol
)

target_link_libraries(statistics_experiments PRIVATE
    helmer_stochastic_sobol
)

# ---------------------------------------------------------------------------
# New cached NUS (Java parent-traversal port)
# ---------------------------------------------------------------------------

set(NEW_CACHED_NUS_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/new-cached-nus/adapter)

add_library(new_cached_nus STATIC
    ${NEW_CACHED_NUS_ADAPTER_DIR}/NewCachedNusPointSet.cpp
)

target_include_directories(new_cached_nus PUBLIC
    ${NEW_CACHED_NUS_ADAPTER_DIR}
)

target_include_directories(new_cached_nus PRIVATE
    ${HELMER_DIR}/sampling
)

target_link_libraries(new_cached_nus PRIVATE
    pbrt_sobol_core
)

add_executable(test_new_cached_nus
    ${NEW_CACHED_NUS_ADAPTER_DIR}/test_adapter.cpp
)

target_include_directories(test_new_cached_nus PRIVATE
    ${HELMER_DIR}/sampling
)

target_link_libraries(test_new_cached_nus PRIVATE
    new_cached_nus
    pbrt_sobol_core
)

add_test(NAME new_cached_nus COMMAND test_new_cached_nus)

target_link_libraries(time_experiments PRIVATE new_cached_nus)
target_link_libraries(statistics_experiments PRIVATE new_cached_nus)

add_executable(test_rqmc_replications tests/test_rqmc_replications.cpp)
target_link_libraries(test_rqmc_replications PRIVATE
    rqmc_common
    helmer_stochastic_sobol
    new_cached_nus
)
add_test(NAME rqmc_replications COMMAND test_rqmc_replications)

# ---------------------------------------------------------------------------
# Cached NUS (sscpp sorted-order port)
# ---------------------------------------------------------------------------

set(CACHED_NUS_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/cached-nus/adapter)

add_library(cached_nus STATIC
    ${CACHED_NUS_ADAPTER_DIR}/CachedNUSPointSet.cpp
)

target_include_directories(cached_nus PUBLIC
    ${CACHED_NUS_ADAPTER_DIR}
)

target_link_libraries(cached_nus PRIVATE
    pbrt_sobol_core
)

# Single-configuration comparison with editable settings and a console report.
add_executable(compare_nus_methods src/experiments/CompareNUSMethods.cpp)
target_link_libraries(compare_nus_methods PRIVATE
    rqmc_common
    cached_nus
    new_cached_nus
    helmer_stochastic_sobol
    friedel_keller_libseq
)

# ---------------------------------------------------------------------------
# Adaptive Regular Tiles / ART-Owen adapter
# ---------------------------------------------------------------------------

set(ADAPTIVE_REGULAR_TILES_ADAPTER_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/methods/adaptive-regular-tiles/adapter)

add_library(adaptive_regular_tiles STATIC
    ${ADAPTIVE_REGULAR_TILES_ADAPTER_DIR}/AdaptiveRegularTiles.cpp
)

target_include_directories(adaptive_regular_tiles PUBLIC
    ${ADAPTIVE_REGULAR_TILES_ADAPTER_DIR}
)

target_link_libraries(adaptive_regular_tiles PUBLIC
    pbrt_sobol_core
)

add_executable(test_adaptive_regular_tiles
    ${ADAPTIVE_REGULAR_TILES_ADAPTER_DIR}/test_adapter.cpp
)

target_link_libraries(test_adaptive_regular_tiles PRIVATE
    adaptive_regular_tiles
)

target_link_libraries(time_experiments PRIVATE
    adaptive_regular_tiles
)

target_link_libraries(statistics_experiments PRIVATE
    adaptive_regular_tiles
)
