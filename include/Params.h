#pragma once

#include <cstdint>
#include <cstddef>

#include "bits_header.h"

namespace xvmt {
namespace details {

template <size_t W>
struct MT19937Params;

template <>
struct MT19937Params<32>
{
    using output_word_t = uint32_t;

    static constexpr BitsGenType s_bitsGenType = BitsGenType::MT32;

    static constexpr size_t s_nBits = 19937;                                                 // Mersenne exponent; period is 2^19937 - 1
    static constexpr size_t s_stateWordBits = 32;                                             // word size w in bits (MT paper notation)
    static constexpr int s_N = s_nBits / s_stateWordBits + (s_nBits % s_stateWordBits != 0);  // 624 - state array length in words
    static constexpr int s_M = 397;                                                          // middle-word offset in the twist recurrence
    static constexpr size_t s_nMatrixBits = s_nBits;                                         // transition matrix dimension in bits

    // tempering constants (MT paper notation: u,d,s,b,t,c,l)
    static constexpr int      s_u = 11;
    static constexpr uint32_t s_d = 0xffffffffUL;              // full mask (no-op AND for 32-bit)
    static constexpr int      s_s = 7;
    static constexpr uint32_t s_b = 0x9d2c5680UL;             // tempering mask b
    static constexpr int      s_t = 15;
    static constexpr uint32_t s_c = 0xefc60000UL;             // tempering mask c
    static constexpr int      s_l = 18;

    // legacy aliases
    static constexpr uint32_t s_temperMask1 = s_b;
    static constexpr uint32_t s_temperMask2 = s_c;

    static constexpr uint32_t s_matrixA = 0x9908b0dfUL;       // companion matrix A lower row (the "a" vector)
    static constexpr uint32_t s_upperMask = 0x80000000UL;     // most-significant (w-r) = 1 bit mask
    static constexpr uint32_t s_lowerMask = 0x7fffffffUL;     // least-significant r = 31 bits mask

    // initialisation constants
    static constexpr int      s_initShift    = 30;
    static constexpr uint32_t s_initMul      = 1812433253UL;
    static constexpr uint32_t s_arrayInitMul1 = 1664525UL;
    static constexpr uint32_t s_arrayInitMul2 = 1566083941UL;
    static constexpr uint32_t s_msb          = 0x80000000UL;
    static constexpr uint32_t s_arrayInitSeed = 19650218UL;

    static constexpr size_t s_n32InOneWord = s_stateWordBits / 32;       // uint32 elements per word (= 1)
    static constexpr size_t s_n32InOneState = s_N * s_n32InOneWord;     // uint32 elements in the full state (= 624)

    static constexpr size_t s_stepOutputWordsLog2 = 0;                  // log2 of output words per generator step
};

template <>
struct MT19937Params<64>
{
    using output_word_t = uint64_t;

    static constexpr BitsGenType s_bitsGenType = BitsGenType::MT64;

    static constexpr size_t s_nBits = 19937;
    static constexpr size_t s_stateWordBits = 64;
    static constexpr int s_N = 312;
    static constexpr int s_M = 156;
    static constexpr size_t s_nMatrixBits = s_nBits;

    // tempering constants
    static constexpr int      s_u = 29;
    static constexpr uint64_t s_d = 0x5555555555555555ULL;
    static constexpr int      s_s = 17;
    static constexpr uint64_t s_b = 0x71D67FFFEDA60000ULL;
    static constexpr int      s_t = 37;
    static constexpr uint64_t s_c = 0xFFF7EEE000000000ULL;
    static constexpr int      s_l = 43;

    static constexpr uint64_t s_matrixA  = 0xB5026F5AA96619E9ULL;
    static constexpr uint64_t s_upperMask = 0xFFFFFFFF80000000ULL;  // upper 33 bits
    static constexpr uint64_t s_lowerMask = 0x000000007FFFFFFFULL;  // lower 31 bits

    // initialisation constants
    static constexpr int      s_initShift    = 62;
    static constexpr uint64_t s_initMul      = 6364136223846793005ULL;
    static constexpr uint64_t s_arrayInitMul1 = 3935559000370003845ULL;
    static constexpr uint64_t s_arrayInitMul2 = 2862933555777941757ULL;
    static constexpr uint64_t s_msb          = 1ULL << 63;
    static constexpr uint64_t s_arrayInitSeed = 19650218ULL;

    static constexpr size_t s_n32InOneWord  = s_stateWordBits / 32;   // 2
    static constexpr size_t s_n32InOneState = s_N * s_n32InOneWord;  // 624

    static constexpr size_t s_stepOutputWordsLog2 = 0;               // log2 of output words per generator step
};

struct SFMT19937Params
{
    static constexpr BitsGenType s_bitsGenType = BitsGenType::SFMT;

    static constexpr size_t s_nBits = 19937;                                                  // Mersenne exponent; period is 2^19937 - 1

    static constexpr size_t s_stateWordBits = 128;                                             // SFMT word size: one 128-bit integer element
    static constexpr int s_N = s_nBits / s_stateWordBits + (s_nBits % s_stateWordBits != 0);   // 156 - state array length in 128-bit words
    static constexpr int s_M = 122;                                                           // POS1 offset in the SFMT recurrence
    static constexpr size_t s_nMatrixBits = s_N * s_stateWordBits;                            // FIXME: confirm this value

    static const size_t s_n32InOneWord = s_stateWordBits / 32;       // uint32 elements per 128-bit word (= 4)
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;     // uint32 elements in the full state (= 624)

    // bit-masks applied to the POS1 element in the SFMT recurrence to ensure maximal period
    static constexpr uint32_t s_SFMT_MSK1 = 0xdfffffefU;
    static constexpr uint32_t s_SFMT_MSK2 = 0xddfecb7fU;
    static constexpr uint32_t s_SFMT_MSK3 = 0xbffaffffU;
    static constexpr uint32_t s_SFMT_MSK4 = 0xbffffff6U;

    static constexpr size_t s_stepOutputWordsLog2 = 2;               // log2 of output words per generator step
};

} // namespace details
} // namespace xvmt
