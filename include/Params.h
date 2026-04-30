#pragma once

#include <cstdint>
#include <cstddef>

namespace xvmt::details {

struct MT19937Params
{
    static constexpr size_t s_nBits = 19937;                                                 // Mersenne exponent; period is 2^19937 - 1
    static constexpr size_t s_wordSizeBits = 32;                                             // word size w in bits (MT paper notation)
    static constexpr int s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);  // 624 - state array length in words
    static constexpr int s_M = 397;                                                          // middle-word offset in the twist recurrence
    static constexpr size_t s_nMatrixBits = s_nBits;                                         // transition matrix dimension in bits

    static constexpr uint32_t s_temperMask1 = 0x9d2c5680UL;   // tempering mask b (MT paper sec.3)
    static constexpr uint32_t s_temperMask2 = 0xefc60000UL;   // tempering mask c (MT paper sec.3)

    static constexpr uint32_t s_matrixA = 0x9908b0dfUL;       // companion matrix A lower row (the "a" vector)
    static constexpr uint32_t s_upperMask = 0x80000000UL;     // most-significant (w-r) = 1 bit mask
    static constexpr uint32_t s_lowerMask = 0x7fffffffUL;     // least-significant r = 31 bits mask

    static constexpr size_t s_n32InOneWord = s_wordSizeBits / 32;       // uint32 elements per word (= 1)
    static constexpr size_t s_n32InOneState = s_N * s_n32InOneWord;     // uint32 elements in the full state (= 624)
};

struct SFMT19937Params
{
    static constexpr size_t s_nBits = 19937;                                                  // Mersenne exponent; period is 2^19937 - 1

    static constexpr size_t s_wordSizeBits = 128;                                             // SFMT word size: one 128-bit integer element
    static constexpr int s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);   // 156 - state array length in 128-bit words
    static constexpr int s_M = 122;                                                           // POS1 offset in the SFMT recurrence
    static constexpr size_t s_nMatrixBits = s_N * s_wordSizeBits;                            // FIXME: confirm this value

    static const size_t s_n32InOneWord = s_wordSizeBits / 32;       // uint32 elements per 128-bit word (= 4)
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;     // uint32 elements in the full state (= 624)

    // bit-masks applied to the POS1 element in the SFMT recurrence to ensure maximal period
    static constexpr uint32_t s_SFMT_MSK1 = 0xdfffffefU;
    static constexpr uint32_t s_SFMT_MSK2 = 0xddfecb7fU;
    static constexpr uint32_t s_SFMT_MSK3 = 0xbffaffffU;
    static constexpr uint32_t s_SFMT_MSK4 = 0xbffffff6U;
};

} // namespace xvmt::details
