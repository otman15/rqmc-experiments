#ifndef ART_H
#define ART_H

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <cmath>
#include <x86intrin.h>

// =============================================================================
// Utils
// =============================================================================

template <typename T>
int findOrInsert(T entry, std::vector<T> &list) {
    int index;
    for (index = 0; index < list.size(); index++) {
        if (list[index] == entry) return index;
    }
    list.push_back(entry);
    return index;
}

// =============================================================================
// Data structures:
// =============================================================================

struct ProductionRule {
    int left, right;
    inline int& operator[](int i) { return ((int *)this)[i]; };
};

typedef std::vector<ProductionRule> Production;

// =============================================================================
// The Thue-Morse word
// =============================================================================

struct {
    inline int getBit (uint32_t index) { return _popcnt32(index) & 1; }
    inline int operator[] (uint32_t index) { return getBit(index); }
    inline uint32_t operator() (uint32_t index, int n) {                        // Retrieves a subword of length n starting at <index>; t[index] is placed in the most-significant bit position.
        uint32_t word = 0;
        for (int i = 0; i < n; i++) {
            word <<= 1;
            word |= getBit(index + i);
        }
        return word;
    }
    inline uint64_t morph(uint32_t word) {                                      // Applies the Thue-Morse morphism to the given word.
        return (
            _pdep_u64( word, 0xAAAAAAAAAAAAAAAAul) |                            // Image of the ones: each one in index i is moved to index 2i
            _pdep_u64(~word, 0x5555555555555555ul)
        );
    }
} t;

// =============================================================================
// ART
// =============================================================================

class ART {
    int idLength;                                                               // Length of identifying words
    int idCount;                                                                // Number of symbols
    uint32_t mask;                                                              // Mask to isolute identifying bits.
    std::vector<uint32_t> idWords;                                              // Actual words for identifying symbols
public:
    ART(int idLength);
    int getIdCount() { return idCount; };
    std::vector<uint32_t> getIdWords() { return idWords; };
    Production getProduction(bool alternate = false);
    char getSymbol(int id);                                                     // An auxiliary function to encode id's consistently as alphabets
};

ART::ART(int idLength) : idLength(idLength) {
    mask = (1 << idLength) - 1;
    int log2idLength2 = std::ceil(std::log2(idLength));                         // Log2 of a power of two no-less that idLength
    int idLength2 = 1 << log2idLength2;                                         // A power of two no less than idLength
    for (int i = 0; i < 7 * idLength2; i++) {                                   // This includes the images of 01, 11, 10, and 00, hence all factors.
        uint32_t idWord = t(i, idLength);
        uint32_t idWordNeg = ~idWord & mask;                                    // Image of idWord with flipped bits
        if (idWord > idWordNeg) std::swap(idWord, idWordNeg);                   // Ensure that idWord starts with 0 and its negative with 1
        findOrInsert(idWord, idWords);                                          // Insert idWord first
        findOrInsert(idWordNeg, idWords);                                       // Insert its negative next
    }
    idCount = idWords.size();
}

Production ART::getProduction(bool alternate) {
    Production production(idCount);
    for (int i = 0; i < idCount; i++) {
        uint32_t idWord = idWords[i];
        if (alternate) idWord = ~idWord;
        uint32_t idLeft  = (t.morph(idWord) >> (idLength    )) & mask;          // Id of leftmost child: morphing doubles the length, then we trim the rightmost half
        uint32_t idRight = (t.morph(idWord) >> (idLength - 1)) & mask;          // Id of rightmost child, starting at the second bit
        int left  = findOrInsert(idLeft, idWords);
        int right = findOrInsert(idRight, idWords);
        if (left >= idCount || right >= idCount) {
            fprintf(stderr, "Id set not closed under morphism\n");
            exit(1);
        }
        production[i] = {left, right};
    }
    return production;
}

char ART::getSymbol(int id) {
    int pairIndex = id / 2;
    char zero = id & 1 ? 'A' : 'a';
    char symbol = pairIndex +  zero;
    return symbol;
}

#endif






