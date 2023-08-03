#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#ifndef SIMD_N_BITS
#   define SIMD_N_BITS 64
#endif

#include <cstdint>
#include <cstddef>

enum VMT19937QueryMode { QM_Scalar, QM_Block16, QM_StateSize };


template <size_t RegisterBitLen = SIMD_N_BITS, VMT19937QueryMode QueryMode = QM_Scalar>
class VMT19937
{
    const static size_t s_regLenBits = RegisterBitLen;
    const static size_t s_regLenWords = s_regLenBits / 32;
    typedef SimdRegister<s_regLenBits> XV;
    typedef SimdRegister<SIMD_N_BITS> XVMax;

    // Period parameters
    static const size_t s_nBits = 19937;
    static const size_t s_N = s_nBits / 32 + (s_nBits % 32 != 0); // 624
    static const size_t s_M = 397;

    static const uint32_t s_rndBlockSize = 64 / sizeof(uint32_t); // exactly one cache line
    
    static const uint32_t s_matrixA = 0x9908b0dfUL;   // constant vector a
    static const uint32_t s_upperMask = 0x80000000UL; // most significant w-r bits
    static const uint32_t s_lowerMask = 0x7fffffffUL; // least significant r bits
    static const uint32_t s_temperMask1 = 0x9d2c5680UL;
    static const uint32_t s_temperMask2 = 0xefc60000UL;

    alignas(64) XV m_state[s_N];    // the array of state vectors

    // This data members is necessary only if QueryMode==QM_Scalar
    alignas(64) uint32_t m_rnd[s_rndBlockSize]; // buffer of tempered numbers with same size as a cache line
    const uint32_t* m_prnd;

    // This data members are redundant if QueryMode==QM_StateSize
    const XV *m_pst, *m_pst_end;    // m_pos==m_pst_end means the state vector has been consumed and need to be regenerated

    static const inline XV v_upperMask = XV(s_upperMask);
    static const inline XV v_lowerMask = XV(s_lowerMask);
    static const inline XV v_matrixA = XV(s_matrixA);
    static const inline XVMax v_temperMask1 = XVMax(s_temperMask1);
    static const inline XVMax v_temperMask2 = XVMax(s_temperMask2);

    template <typename U>
    static FORCE_INLINE U temper(U y, U mask1, U mask2)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & mask1);
        y = y ^ ((y << 15) & mask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <bool Aligned>
    static FORCE_INLINE void temperRefillBlock(const XV *&st, uint32_t *dst)
    {
        typedef SimdRegister<SIMD_N_BITS> XVmax;

        const XVmax mask1(v_temperMask1);
        const XVmax mask2(v_temperMask2);
        XVMax* prnd = (XVMax*)dst;
        const XVMax* pst = (XVMax*)st;
        for (size_t i = 0; i < s_rndBlockSize * sizeof(uint32_t) / sizeof(XVmax); ++i) {
            XVmax tmp = temper(*pst++, mask1, mask2);
            tmp.template store<Aligned>((uint32_t*)(prnd++));
        }
        st = (const XV*)(pst);
    }

    static FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm, const XV& upperMask, const XV& lowerMask, const XV& matrixA)
    {
        XV y = (s & upperMask) | (sp & lowerMask);
        // y and sp are either both even or both odd,
        // hence in the next line we can check if sp is odd
        // so that the operation is independent on the clauclation of y
        // and the compiler is free to rearrange the code
        XV r = sm ^ (y >> 1) ^ sp.ifOddValueElseZero(matrixA);
        return r;
    }

    void FORCE_INLINE refill()
    {
        XV* stCur = m_state;
        XV* stNxt = m_state + 1;
        const XV* stMid = m_state + s_M;
        const XV* stEnd = m_state + s_N;

        // Load the variables in registers and passes them to the function as argument,
        // to avoid to re-read the static variables from memory at every iteration
        // Note that since the function is forced inline, the function arguments
        // do no need to be loaded on the stack.
        const XV upperMask(v_upperMask);
        const XV lowerMask(v_lowerMask);
        const XV matrixA(v_matrixA);

        //size_t kk;
        XV cur = *stCur;
        for (; stMid != stEnd; ++stMid, stCur = stNxt++) {
            XV nextp = *stNxt;
            *stCur = advance1(cur, nextp, *stMid, upperMask, lowerMask, matrixA);
            cur = nextp;
        }
        for (stMid = m_state; stNxt != stEnd; ++stMid, stCur = stNxt++) {
            XV nextp = *stNxt;
            *stCur = advance1(cur, nextp, *stMid, upperMask, lowerMask, matrixA);
            cur = nextp;
        }
        {
            *stCur = advance1(cur, m_state[0], *stMid, upperMask, lowerMask, matrixA);
        }

        m_pst = m_state;
    }

    uint32_t& scalarState(uint32_t scalarIndex)
    {
        return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    }

    uint32_t scalarState(uint32_t scalarIndex) const
    {
        return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    }

    // extract one of the interleaved state vectors, shift it left by 31 bits and save it to dst
    void stateToVector(size_t stateIndex, uint32_t *pdst) const
    {
        const uint32_t* pstate = (const uint32_t*)m_state;
        pdst[0] = pstate[stateIndex] >> 31;
        for (size_t i = 1; i < s_N; ++i) {
            uint32_t word = pstate[i * s_regLenWords + stateIndex];
            pdst[i - 1] |= word << 1;
            pdst[i] = word >> 31;
        }
    }

    // shift vector psr to the right by 31 bit and store into the interleaved elements of the state vector
    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        uint32_t* pstate = (uint32_t*)m_state;
        const uint32_t* pw = (const uint32_t*)psrc;
        pstate[stateIndex] = 0;
        size_t w;
        for (w = 0; w < s_N - 1; ++w) {
            uint32_t word = pw[w];
            pstate[w * s_regLenWords + stateIndex] |= word << 31;
            pstate[(w + 1) * s_regLenWords + stateIndex] = word >> 1;
        }
        pstate[w * s_regLenWords + stateIndex] |= pw[w] << 31;
    }

    void fillOtherStates(size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
    {
        uint32_t* pstate = (uint32_t*)m_state;

        // temporary workspace matrix
        BinaryMatrix<2, s_nBits> tmp;

        if (commonJumpRepeat) {
            MYASSERT(commonJump, "commonJump is required when commnJumpRepeat>0");

            // copy state to the first row shifting all bits to the left by 31
            stateToVector(0, (uint32_t*)tmp.rowBegin(0));

            for (size_t i = 0; i < commonJumpRepeat; ++i)
                commonJump->multiplyByColumn(tmp.rowBegin((i+1) % 2), tmp.rowBegin(i % 2));

            // copy to the state vector shifting all bits to the right by 31
            vectorToState(0, (const uint32_t*) tmp.rowBegin(commonJumpRepeat % 2));
        }

        if (s_regLenWords > 1) {
            if (sequentialJump) {

                // perform jump ahead of the s_regLenWords states
                // State_0 = State_0
                // State_1 = Jump x State_0
                // State_2 = Jump x State_1
                // ...

                // copy state to the first row shifting all bits to the left by 31
                stateToVector(0, (uint32_t*)tmp.rowBegin(0));

                for (size_t s = 1; s < s_regLenWords; ++s) {
                    // multiply all rows by state s and store the result in pres

                    const uint8_t* psrc = (uint8_t*)tmp.rowBegin((s + 1) % 2);
                    uint8_t* pdst = (uint8_t*)tmp.rowBegin(s % 2);

                    sequentialJump->multiplyByColumn(pdst, psrc);

                    // copy to the state vector shifting all bits to the right by 31
                    vectorToState(s, (const uint32_t*)pdst);
                }
            }
            else {
                for (size_t w = 0; w < s_N; ++w)
                    for (size_t j = 1; j < s_regLenWords; ++j)
                        pstate[w * s_regLenWords + j] = pstate[w * s_regLenWords];
            }
        }

        m_pst = m_pst_end;
        m_prnd = (const uint32_t *)(((uint8_t *) m_rnd) + sizeof(m_rnd));
    }

    // initializes the first state with a seed
    void __reinit(uint32_t s)
    {
        const uint32_t mask = uint32_t(1812433253UL);
        uint32_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < s_N; i++)
            prev = scalarState(i) = (mask * (prev ^ (prev >> 30)) + i);
    }

public:

    const static size_t s_qryBlkSize = (QueryMode == QM_Scalar) ? 1 : (QueryMode == QM_Block16) ? 16 : s_N * s_regLenWords;

    // constructors
    VMT19937()
        : m_pst(nullptr)
        , m_pst_end(m_state+s_N)
        , m_prnd(nullptr)
    {}

    VMT19937(uint32_t seed, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
        : VMT19937()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    VMT19937(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
        : VMT19937()
    {
        reinit(seeds, n_seeds, commonJumpRepeat, commonJump, sequentialJump);
    }

    // initializes m_state[s_N] with a seed
    void reinit(uint32_t s, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
    {
        __reinit(s);
        fillOtherStates(commonJumpRepeat, commonJump, sequentialJump);
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    void reinit(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
    {
        __reinit(uint32_t(19650218));
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

        fillOtherStates(commonJumpRepeat, commonJump, sequentialJump);
    }


    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        static_assert(QueryMode == QM_Scalar);

        if ((reinterpret_cast<intptr_t>(m_prnd) % sizeof(m_rnd)) != 0)
            return *m_prnd++;

        if (m_pst != m_pst_end)
            /* do nothing*/; // most likely case first
        else
            refill();

        temperRefillBlock<true>(m_pst, m_rnd);
        m_prnd = m_rnd;

        return *m_prnd++;
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QueryMode == QM_Block16);

        if (m_pst != m_pst_end)
            /* do nothing*/; // most likely case first
        else
            refill();
        temperRefillBlock<false>(m_pst, dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        static_assert(QueryMode == QM_StateSize);

        refill();
        const XV* pst = m_state;
        for (size_t i = 0; i < sizeof(m_state) / sizeof(m_rnd); ++i, dst += s_rndBlockSize)
            temperRefillBlock<false>(pst, dst);
    }

    // generates a random number on [0,0x7fffffff]-interval
    static uint32_t convert_uint31(uint32_t rnd)
    {
        return (uint32_t)(rnd >> 1);
    }

    // generates a random number on [0,1]-real-interval
    static double convert_real1(uint32_t rnd)
    {
        return rnd * (1.0 / 4294967295.0);
        // divided by 2^32-1 
    }

    // generates a random number on [0,1)-real interval 
    static double convert_real2(uint32_t rnd)
    {
        return rnd * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on (0,1)-real-interval
    static double convert_real3(uint32_t rnd)
    {
        return (((double)rnd) + 0.5) * (1.0 / 4294967296.0);
        // divided by 2^32 
    }

    // generates a random number on [0,1) with 53-bit resolution
    static double convert_res53(uint32_t rnd1, uint32_t rnd2)
    {
        uint32_t a = rnd1 >> 5, b = rnd2 >> 6;
        return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }
};
