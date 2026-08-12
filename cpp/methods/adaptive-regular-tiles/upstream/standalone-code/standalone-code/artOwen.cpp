#include "art.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <random>
#include <stdint.h>
#include <getopt.h>


// =============================================================================
// Data tyoes
// =============================================================================

struct Point {
    uint32_t x, y;
};

// typedef uint32_t unsigned;
typedef Production Prod;
const int BITS = 32;

typedef std::vector<Point> Points;

// =============================================================================
// Utils:
// =============================================================================

std::random_device rd;                                                          // Will be used to obtain a seed for the random number engine
std::mt19937 rnd(rd());

inline void shuffle(std::vector<int> &list) {                                   // Generate a permutation of {0..N-1}
    uint32_t N = list.size();
    for (uint32_t i = 0; i < N; i++) list[i] = i;
    for (uint32_t i = N - 1; i > 0; i--) {                                      // Iterate down the list, and swap each place with another place down, inclusive of same place.
        uint32_t r = rnd() % (i + 1);                                           // (i + 1) means that the same place is included.
        std::swap(list[i], list[r]);
    }
}

// =============================================================================
// Sobol
// =============================================================================

inline uint32_t bitReverse(uint32_t i) {
    i = (i << 16) | (i >> 16);
    i = ((i & 0x00FF00FF) << 8) | ((i & 0xFF00FF00) >> 8);
    i = ((i & 0x0F0F0F0F) << 4) | ((i & 0xF0F0F0F0) >> 4);
    i = ((i & 0x33333333) << 2) | ((i & 0xCCCCCCCC) >> 2);
    i = ((i & 0x55555555) << 1) | ((i & 0xAAAAAAAA) >> 1);
    return i;
}

inline uint32_t sobol(uint32_t r) {                                             // Input is a bit-reversed uint32_t
    r ^= (r & 0x0000FFFF) << 16;
    r ^= (r & 0x00FF00FF) <<  8;
    r ^= (r & 0x0F0F0F0F) <<  4;
    r ^= (r & 0x33333333) <<  2;
    r ^= (r & 0x55555555) <<  1;
    return r;
}

// =============================================================================
// ART-Owen scrambling
// =============================================================================

inline unsigned scrmbl(unsigned xIn, Prod &prod, unsigned *data, int id) {
    unsigned xOut = xIn;
    for (int i = 0; i < BITS; i++) {
        xOut ^=  data[id] >> i;
        id = prod[id][(xIn >> ((BITS - 1) - i)) & 1];
    }
    return xOut;
}

inline unsigned unscrmbl(unsigned xOut, Prod &prod, unsigned *data, int id) {
    unsigned xIn = xOut;
    for (int i = 0; i < BITS; i++) {
        xIn ^=  data[id] >> i;
        id = prod[id][(xIn >> ((BITS - 1) - i)) & 1];
    }
    return xIn;
}

inline unsigned scrmbl01(unsigned xIn, unsigned *data, int id = 0) {
    unsigned xOut = xIn;
    for (int i = 0; i < BITS; i++) {
        xOut ^=  data[id] >> i;
        id ^= (xIn >> ((BITS - 1) - i)) & 1;
    }
    return xOut;
}


// =============================================================================
// Generate ordered grammar
// =============================================================================

Production ordered(int alphabetSize) {
    Production prod(alphabetSize, {-1, -1});
    for (int i = 2; i <= alphabetSize; i++) {                                   // The trick is to number the nodes from 1
        int parent = i >> 1;
        prod[parent - 1][i & 1] = i - 1;
    }
    for (int i = 0; i < alphabetSize; i++) {
        if (prod[i].left  < 0) prod[i].left  = rnd() % alphabetSize;
        if (prod[i].right < 0) prod[i].right = rnd() % alphabetSize;
    }
    return prod;
}

// =============================================================================


const char *USAGE_MESSAGE = "Usage: %s [options] [outputfile]\n"
"Options:\n"
"   -q                      Suppress messages\n"
"   -m <log2N>              Size of generated sets, default is 12 (4k)\n"
"   -a <idLength>           Default is 2\n"
"   -R                      Generate a random grammar instead of ART\n"
"   -A                      Morph alternately between t and its complement\n"
"   -e                      Print production as eps array\n"
"   -g                      Read grammar from stdin\n"
"   -o                      Generate ordered grammar\n"
;

int main(int argc,char **argv) {
    int opt;
    bool verbose = true;
    int idLength = 2;
    bool useRandomGrammar = false;
    int m = 12;
    bool alternate = false;
    bool epsProduction = false;
    bool readGrammar = false;
    bool useOrderedGrammar = false;
    int rootId = 0;
    while ((opt = getopt(argc, argv, "qa:Rm:Aego")) != -1) {
        switch (opt) {
            case 'q': verbose = false; break;
            case 'a': idLength = atoi(optarg); break;
            case 'R': useRandomGrammar = true; break;
            case 'm': m = atoi(optarg); break;
            case 'A': alternate = true; break;
            case 'e': epsProduction = true; break;
            case 'g': readGrammar = true; break;
            case 'o': useOrderedGrammar = true; break;
            default: fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
        }
    }
    if (optind > argc) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]); exit(1);
    }
    unsigned N = 1 << m;
    Points p(N);

    ART art(idLength);
    int alphabetSize = art.getIdCount();
    std::vector<unsigned> scrambleCodeX(alphabetSize), scrambleCodeY(alphabetSize);
    Production production = art.getProduction(alternate);
    // ============================
    // Random grammar
    // ============================
    if (useRandomGrammar) {
        for (int id = 0; id < alphabetSize; id++) {
            production[id] = {rnd() % alphabetSize, rnd() % alphabetSize};
        }
    }
    // ============================
    // Ordered grammar
    // ============================
    if (useOrderedGrammar) {
        production = ordered(alphabetSize);
    }
    // ============================
    // Read grammar from stdin:
    // ============================
    if (readGrammar) {
        fprintf(stderr, "Please enter rules seprated by white space:\n");
        for (int id = 0; id < alphabetSize; id++) {
            int left, right;
            fprintf(stderr, "%2d --> \n");
            scanf(" %d %d", &left, &right);
            production[id] = {left, right};
        }
    }
    // ============================
    // ============================
    if (verbose) {
        fprintf(stderr, "Alphabet Size = %d\n", alphabetSize);
        for (int i = 0; i < alphabetSize; i++) {
            fprintf(
                stderr, "%c --> %c%c\n",
                art.getSymbol(i),
                art.getSymbol(production[i].left),
                art.getSymbol(production[i].right)
            );
        }
    }
    if (epsProduction) {
        printf("/production [");
        for (int i = 0; i < alphabetSize; i++) {
            printf("[%d %d]", production[i].left, production[i].right);
        }
        printf("] def\n");
    }
    // =========================================================================
    // Save point sets
    // =========================================================================
    if (optind < argc) {                                                        // File names supplied?
        int n = argc - optind;
        double clocks = 0;
        for (int fileNo = 0; fileNo < n; fileNo++) {
            clock_t t0 = clock();
            for (int i = 0; i < alphabetSize; i++) {
                scrambleCodeX[i] = rnd();
                scrambleCodeY[i] = rnd();
            }
            if (alphabetSize == 2) {
                for (unsigned i = 0; i < N; i++) {                                          // Iterate through bit-reversed indices
                    unsigned x = bitReverse(i);
                    unsigned y = sobol(x);
                    p[i] = {
                        scrmbl01(x, scrambleCodeX.data(), rootId),
                        scrmbl01(y, scrambleCodeY.data(), rootId)
                    };
                }
            }
            else {
                for (unsigned i = 0; i < N; i++) {                                          // Iterate through bit-reversed indices
                    unsigned x = bitReverse(i);
                    unsigned y = sobol(x);
                    p[i] = {
                        scrmbl(x, production, scrambleCodeX.data(), rootId),
                        scrmbl(y, production, scrambleCodeY.data(), rootId)
                    };
                }
            }
            clocks += clock() - t0;
            FILE *file = fopen(argv[optind + fileNo], "w");
            const double inv = 1.0 / (1ull << BITS);
            fprintf(file, "%d\n", N);
            for (int i = 0; i < N; i++) {
                fprintf(file, "%0.18f %0.18f\n", inv * p[i].x, inv * p[i].y);
            }
            fclose(file);
        }
        if (verbose) {
            fprintf(
                stderr, "Average generation rate is %4dM points per second)\n",
                int(n * N * double(CLOCKS_PER_SEC) / (clocks * 1000000.0))
            );
        }
    }
}





