#pragma once

/* 
   Converted to C++ from the original code downloaded from:
   mt19937ar.c
   http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html
   and vectorized using SIMD.
*/

#include "SIMD.h"

#ifndef SIMD_N_BITS
#   define SIMD_N_BITS 64
#endif

#include <cstdint>
#include <cstddef>

template <size_t BitsLen = SIMD_N_BITS>
class MT19937SIMD
{
    const static size_t s_regLenBits = BitsLen;
    const static size_t s_regLenWords = s_regLenBits / 32;
    typedef SimdRegister<s_regLenBits> XV;

    // Period parameters
    static const size_t s_N = 624;
    static const size_t s_M = 397;
    static const uint32_t s_matrixA = 0x9908b0dfUL;   // constant vector a
    static const int32_t s_upperMask = 0x80000000UL; // most significant w-r bits
    static const int32_t s_lowerMask = 0x7fffffffUL; // least significant r bits
    static const uint32_t s_temperMask1 = 0x9d2c5680UL;
    static const uint32_t s_temperMask2 = 0xefc60000UL;

    XV m_state[s_N];  // the array os state vectors
    XV m_u32[s_N]; // a cache of uniform discrete random numbers in the range [0,0xffffffff]
    uint32_t m_pos = s_N + 1;    // m_pos==s_N+1 means m_state[s_N] is not initialized

    struct Loop
    {
        const XV v_matrixA;
        const XV v_upperMask;
        const XV v_lowerMask;
        const XV v_temperMask1;
        const XV v_temperMask2;
        const XV v_one;
        const XV v_zero;

        Loop()
            : v_matrixA(s_matrixA)
            , v_upperMask(s_upperMask)
            , v_lowerMask(s_lowerMask)
            , v_temperMask1(s_temperMask1)
            , v_temperMask2(s_temperMask2)
            , v_one(uint32_t(0x1))
            , v_zero(XV::zero())
        {
        }

        FORCE_INLINE XV temper(XV y)
        {
            y = y ^ (y >> 11);
            y = y ^ ((y << 7) & v_temperMask1);
            y = y ^ ((y << 15) & v_temperMask2);
            y = y ^ (y >> 18);
            return y;
        }

        FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm)
        {
            XV y = (s & v_upperMask) | (sp & v_lowerMask);
            XV r = sm ^ (y >> 1) ^ (((y & v_one) > v_zero) & v_matrixA);
            return r;
        }
    };

    void refill()
    {

        Loop loop;

        size_t kk;
        XV cur = m_state[0];
        for (kk = 0; kk < s_N - s_M; kk++) {
            XV nextp = m_state[kk + 1];
            XV tmp = loop.advance1(cur, nextp, m_state[kk + s_M]);;
            m_state[kk] = tmp;
            m_u32[kk] = loop.temper(tmp);
            cur = nextp;
        }
        for (; kk < s_N - 1; kk++) {
            XV nextp = m_state[kk + 1];
            XV tmp = loop.advance1(cur, nextp, m_state[kk - (s_N - s_M)]);
            m_state[kk] = tmp;
            m_u32[kk] = loop.temper(tmp);
            cur = nextp;
        }
        {
            XV tmp = loop.advance1(cur, m_state[0], m_state[s_M - 1]);
            m_state[s_N - 1] = tmp;
            m_u32[s_N - 1] = loop.temper(tmp);
        }

        m_pos = 0;
    }

    uint32_t& scalarState(uint32_t scalarIndex)
    {
        return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    }

public:
    // constructors
    MT19937SIMD() : m_pos(s_N + 1) {}

    MT19937SIMD(uint32_t seed)
    {
        reinit(seed);
    }

    MT19937SIMD(const uint32_t seeds[], uint32_t n_seeds)
    {
        reinit(seeds, n_seeds);
    }

    // initializes m_state[s_N] with a seed
    void reinit(uint32_t s, bool fillOthers = true)
    {
        const uint32_t mask = uint32_t(1812433253UL);
        uint32_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < s_N; i++)
            prev = scalarState(i) = (mask * (prev ^ (prev >> 30)) + i);
        m_pos = s_N * s_regLenWords;
        if(fillOthers)
            fillOtherStates();
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    // slight change for C++, 2004/2/26
    void reinit(const uint32_t seeds[], uint32_t n_seeds)
    {
        reinit(uint32_t(19650218), false);
        uint32_t i = 1, j = 0;
        uint32_t k = (s_N > n_seeds ? s_N : n_seeds);
        for (; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i-1) ^ (scalarState(i-1) >> 30)) * uint32_t(1664525)))
                + seeds[j] + j; // non linear
            //m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++; j++;
            if (i >= s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
            if (j >= n_seeds) j = 0;
        }
        for (k = s_N - 1; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * uint32_t(1566083941)))
                - i; // non linear 
            //m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++;
            if (i >= s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
        }

        scalarState(0) = uint32_t(0x80000000); // MSB is 1; assuring non-zero initial array

        fillOtherStates();
    }

    void fillOtherStates()
    {
        uint32_t* p = (uint32_t*)m_state;
        for (size_t i = 0; i < s_N; ++i)
            for (size_t j = 1; j < s_regLenWords; ++j)
                p[i * s_regLenWords + j] = p[i * s_regLenWords];
    }

    // generates a random number on [0,0xffffffff] interval
    uint32_t genrand_uint32()
    {
        if (m_pos < s_N * s_regLenWords)
            /* do nothing*/;
        else
            refill();
        return ((const uint32_t *) m_u32)[m_pos++];
    }

    // generates a random number on [0,0x7fffffff]-interval
    uint32_t genrand_uint31(void)
    {
        return (long)(genrand_uint32() >> 1);
    }

    // generates a random number on [0,1]-real-interval
    double genrand_real1(void)
    {
        return genrand_uint32() * (1.0 / 4294967295.0);
        // divided by 2^32-1 
    }

    // generates a random number on [0,1)-real interval 
    double genrand_real2(void)
    {
        return genrand_uint32() * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on (0,1)-real-interval
    double genrand_real3(void)
    {
        return (((double)genrand_uint32()) + 0.5) * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on [0,1) with 53-bit resolution
    double genrand_res53(void)
    {
        uint32_t a = genrand_uint32() >> 5, b = genrand_uint32() >> 6;
        return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }
    // These real versions are due to Isaku Wada, 2002/01/09 added
};
