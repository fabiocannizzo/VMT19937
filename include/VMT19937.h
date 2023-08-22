#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>

#ifndef VMT19937_STATIC_CONST
#   define VMT19937_STATIC_CONST 0
#endif

namespace Details {

template <size_t RegisterBitLen = SIMD_N_BITS>
class VMT19937Base
{

    static const size_t s_nBits = 19937;
    static const size_t s_wordSizeBits = 32;

    static_assert(RegisterBitLen >= s_wordSizeBits);

    static const size_t s_n32inReg = RegisterBitLen / 32;

    static const int s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);     // 624
    static_assert(s_N == 624);

public:
    static const size_t s_regLenBits = RegisterBitLen;
    static const size_t s_nStates = RegisterBitLen / s_wordSizeBits;
    static const size_t s_n32InOneWord = s_wordSizeBits / 32;            // 4
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;          // 624
    const static size_t s_n32InFullState = s_n32InOneState * s_nStates;  // 624 * nStates

    typedef MT19937Matrix matrix_t;

private:
    const static size_t s_regLenWords = s_regLenBits / s_wordSizeBits;  // FIXME: review this definition
    typedef SimdRegister<s_regLenBits> XV;
    typedef SimdRegister<SIMD_N_BITS> XVMax;

    // Period parameters
    static const size_t s_M = 397;

    static const uint32_t s_rndBlockSize = 64 / sizeof(uint32_t); // exactly one cache line

    alignas(64) XV m_state[s_N];    // the array of state vectors

    // This data members is necessary only if QueryMode==QM_Scalar
    alignas(64) uint32_t m_rnd[s_rndBlockSize]; // buffer of tempered numbers with same size as a cache line
    const uint32_t* m_prnd;

    // This data members are redundant if QueryMode==QM_StateSize
    const XV *m_pst, *m_pst_end;    // m_pos==m_pst_end means the state vector has been consumed and need to be regenerated

private:

    template <typename XVI>
    class TemperCst
    {
        static const uint32_t s_temperMask1 = 0x9d2c5680UL;
        static const uint32_t s_temperMask2 = 0xefc60000UL;
    public:
        TemperCst() : m_mask1(s_temperMask1), m_mask2(s_temperMask2) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    class RefillCst
    {
        static const uint32_t s_matrixA = 0x9908b0dfUL;   // constant vector a
        static const uint32_t s_upperMask = 0x80000000UL; // most significant w-r bits
        static const uint32_t s_lowerMask = 0x7fffffffUL; // least significant r bits
    public:
        RefillCst() : m_upperMask(s_upperMask), m_lowerMask(s_lowerMask), m_matrixA(s_matrixA) {}
        const XV m_upperMask;
        const XV m_lowerMask;
        const XV m_matrixA;
    };
#if (VMT19937_STATIC_CONST==1)
    const static RefillCst s_refillCst;
#endif

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

        TemperCst<XVmax> cst{};
        XVmax* prnd = (XVmax*)dst;
        const XVMax* pst = (XVMax*)st;
        for (size_t i = 0; i < s_rndBlockSize * sizeof(uint32_t) / sizeof(XVmax); ++i) {
            XVmax tmp = temper(*pst++, cst.m_mask1, cst.m_mask2);
            tmp.template store<Aligned>((uint32_t*)(prnd++));
        }
        st = (const XV*)(pst);
    }

    static FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm, const RefillCst& masks)
    {
        XV y = (s & masks.m_upperMask) | (sp & masks.m_lowerMask);
        // y and sp are either both even or both odd,
        // hence in the next line we can check if sp is odd
        // so that the operation is independent on the calculation of y
        // and the compiler is free to rearrange the code
        XV r = sm ^ (y >> 1) ^ sp.ifOddValueElseZero(masks.m_matrixA);
        return r;
    }

    template <int nIter, int J0, int J1, int JM>
    static FORCE_INLINE void unroll(XV* p, XV& pJ0, const RefillCst& masks)
    {
        if constexpr (nIter > 0) {
#ifdef DEBUG
            if (!pJ0.eq(p[J0]))
                THROW("how did this happen?");
#endif
            XV pJ1(p[J1]);
            p[J0] = advance1(pJ0, pJ1, p[JM], masks);
            pJ0 = pJ1;
            unroll<nIter - 1, J0 + 1, J1 + 1, JM + 1>(p, pJ0, masks);
        }
    }

    template <int nUnroll, int nIter, int J1, int JM>
    static FORCE_INLINE void advanceLoop(XV*& p, XV& pJ0, const RefillCst& masks)
    {
        if constexpr (nIter >= nUnroll) {
            auto pend = p + (nIter / nUnroll) * nUnroll;
            // unroll the loop in blocks of nUnroll
            do {
                unroll<nUnroll, 0, J1, JM>(p, pJ0, masks);
                p += nUnroll;
            } while (p != pend);
        }
        const size_t nResIter = nIter % nUnroll;
        if constexpr (nResIter) {
            unroll<nResIter, 0, J1, JM>(p, pJ0, masks);
            p += nResIter;
        }
    }

    void NO_INLINE refill()
    {
        XV* stCur = m_state;

        static const int N = s_N;
        static const int M = s_M;
        static_assert(N == 624 && M == 397, "unrolling designed for these parameters");

        // Create local copy of the constants and pass them to the function as arguments.
        // Since all functions invoked from here are forced inline, the function arguments
        // will not be passed as arguments via the stack, but reside in CPU registers
#if (VMT19937_STATIC_CONST==1)
        const RefillCst masks{ s_refillCst };  // use copy constructor
#else
        const RefillCst masks;  // use default constructor
#endif
        XV pJ0 = stCur[0];

        // unroll first part of the loop (N-M) iterations
        advanceLoop<4, N - M, 1, M>(stCur, pJ0, masks);

        // unroll second part of the loop (M-1) iterations
        advanceLoop<4, M - 1, 1, M - N>(stCur, pJ0, masks);

        // last iteration
        advanceLoop<1, 1, 1 - N, M - N>(stCur, pJ0, masks);

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

        if constexpr (s_regLenWords > 1) {
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

    // constructors
    VMT19937Base()
        : m_prnd(nullptr)
        , m_pst(nullptr)
        , m_pst_end(m_state+s_N)
    {}

    VMT19937Base(uint32_t seed, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
        : VMT19937Base()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    VMT19937Base(const uint32_t *seeds, uint32_t nSeeds, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
        : VMT19937Base()
    {
        reinit(seeds, nSeeds, commonJumpRepeat, commonJump, sequentialJump);
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
    void reinit(const uint32_t *seeds, uint32_t nSeeds, size_t commonJumpRepeat, const BinaryMatrix<s_nBits>* commonJump, const BinaryMatrix<s_nBits>* sequentialJump)
    {
        __reinit(uint32_t(19650218));
        uint32_t i = 1, j = 0;
        uint32_t k = (s_N > nSeeds ? s_N : nSeeds);
        for (; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i-1) ^ (scalarState(i-1) >> 30)) * uint32_t(1664525)))
                + seeds[j] + j; // non linear
            //m_state[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
            i++; j++;
            if (i >= s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
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
        refill();
        const XV* pst = m_state;
        for (size_t i = 0; i < s_n32InFullState / s_rndBlockSize; ++i, dst += s_rndBlockSize)
            temperRefillBlock<false>(pst, dst);
    }
};

#if (VMT19937_STATIC_CONST==1)
template <size_t RegisterBitLen, VRandGenQueryMode QueryMode>
const typename VMT19937<RegisterBitLen, QueryMode>::RefillCst VMT19937<RegisterBitLen, QueryMode>::s_refillCst;
#endif

} // namespace Details
