#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>

#ifndef VMT19937_STATIC_CONST
#   define VMT19937_STATIC_CONST 0
#endif

namespace Details {

template <size_t RegisterBitLen, size_t RegisterBitLenImpl>
class VMT19937Base
{

    static const size_t s_nBits = MT19937Params::s_nBits;
    static const size_t s_wordSizeBits = MT19937Params::s_wordSizeBits;

    static_assert(RegisterBitLen >= s_wordSizeBits);

public:

    static const int s_N = MT19937Params::s_N;     // 624
    static_assert(s_N == 624);

    static const size_t s_regLenBits = RegisterBitLen;
    static const size_t s_regLenImplBits = RegisterBitLenImpl;
    static const size_t s_nStates = RegisterBitLen / s_wordSizeBits;
    static const size_t s_n32inReg = RegisterBitLen / 32;
    static const size_t s_n32InOneWord = s_wordSizeBits / 32;            // 1
    static const size_t s_n32InOneState = s_N * s_n32InOneWord;          // 624
    const static size_t s_n32InFullState = s_n32InOneState * s_nStates;  // 624 * nStates
    const static size_t s_nMatrixBits = MT19937Params::s_nMatrixBits;

    typedef MT19937Matrix matrix_t;

private:
    const static size_t s_regLenWords = s_regLenBits / s_wordSizeBits;  // FIXME: review this definition
    typedef SimdRegister<s_regLenBits, RegisterBitLenImpl> XV;

    static const uint32_t s_n32InRndCache = 64 / sizeof(uint32_t); // exactly one cache line

protected:
    alignas(64) uint32_t m_state[s_N * s_n32inReg];    // the array of state vectors

private:
    // This data members is necessary only if QueryMode==QM_Scalar
    alignas(64) uint32_t m_rnd[s_n32InRndCache]; // buffer of tempered numbers with same size as a cache line
    const uint32_t* m_prnd, * const m_rndEnd;

    // This data members are redundant if QueryMode==QM_StateSize
    const uint32_t*m_pst, * const m_pstEnd;    // m_pos==m_pstEnd means the state vector has been consumed and need to be regenerated

private:

    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask1(MT19937Params::s_temperMask1), m_mask2(MT19937Params::s_temperMask2) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    struct RefillCst
    {
        RefillCst() : m_upperMask(MT19937Params::s_upperMask), m_lowerMask(MT19937Params::s_lowerMask), m_matrixA(MT19937Params::s_matrixA) {}
        const XV m_upperMask;
        const XV m_lowerMask;
        const XV m_matrixA;
    };
#if (VMT19937_STATIC_CONST==1)
    const static RefillCst s_refillCst;
#endif

    template <typename XVI>
    static FORCE_INLINE XVI temper(XVI y, const TemperCst<XVI>& masks)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & masks.m_mask1);
        y = y ^ ((y << 15) & masks.m_mask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <bool Aligned>
    static FORCE_INLINE void temperRefillBlock(const uint32_t *&st_, uint32_t *dst)
    {
        typedef SimdRegister<SIMD_N_BITS, SIMD_N_BITS> XVmax;
        const size_t n32PerIteration = sizeof(XVmax) / sizeof(uint32_t);

        TemperCst<XVmax> cst{};
        const uint32_t* st = st_;
        for (size_t i = 0; i < s_n32InRndCache * sizeof(uint32_t) / sizeof(XVmax); ++i) {
            XVmax tmp = temper(XVmax(st), cst);
            tmp.template store<Aligned>((uint32_t*)(dst));
            dst += n32PerIteration;
            st += n32PerIteration;
        }
        st_ = st;
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
    static FORCE_INLINE void unroll(uint32_t* p, XV& x0, const RefillCst& masks)
    {
        if constexpr (nIter > 0) {
#ifdef DEBUG
            if (!pJ0.eq(p[J0]))
                THROW("how did this happen?");
#endif
            XV x1(p + J1 * s_n32inReg);
            XV xM(p + JM * s_n32inReg);
            XV tmp = advance1(x0, x1, xM, masks);
            tmp.template store<true>(p + J0 * s_n32inReg);
            x0 = x1;
            unroll<nIter - 1, J0 + 1, J1 + 1, JM + 1>(p, x0, masks);
        }
    }

    template <int nUnroll, int nIter, int J1, int JM>
    static FORCE_INLINE void advanceLoop(uint32_t*& p, XV& x0, const RefillCst& masks)
    {
         const size_t nBlockIter = nIter / nUnroll;
         const size_t nResIter = nIter % nUnroll;
        if constexpr (nBlockIter > 0) {
            auto pend = p + nBlockIter * nUnroll * s_n32inReg;
            // unroll the loop in blocks of nUnroll
            do {
                unroll<nUnroll, 0, J1, JM>(p, x0, masks);
                p += nUnroll * s_n32inReg;
            } while (p != pend);
        }
        if constexpr (nResIter) {
            unroll<nResIter, 0, J1, JM>(p, x0, masks);
            p += nResIter * s_n32inReg;
        }
    }

    void NO_INLINE refill()
    {
        uint32_t* stCur = m_state;

        const int N = s_N;
        const int M = Details::MT19937Params::s_M;
        static_assert(N == 624 && M == 397, "unrolling designed for these parameters");

        // Create local copy of the constants and pass them to the function as arguments.
        // Since all functions invoked from here are forced inline, the function arguments
        // will not be passed as arguments via the stack, but reside in CPU registers
#if (VMT19937_STATIC_CONST==1)
        const RefillCst masks{ s_refillCst };  // use copy constructor
#else
        const RefillCst masks;  // use default constructor
#endif
        XV x0(stCur);

        // unroll first part of the loop (N-M) iterations
        advanceLoop<4, N - M, 1, M>(stCur, x0, masks);

        // unroll second part of the loop (M-1) iterations
        advanceLoop<4, M - 1, 1, M - N>(stCur, x0, masks);

        // last iteration
        advanceLoop<1, 1, 1 - N, M - N>(stCur, x0, masks);

        m_pst = begin();
    }

    const uint32_t* begin() const
    {
        return m_state;
    }

    const uint32_t* end() const
    {
        return m_pstEnd;
    }

    const uint32_t* beginRnd() const
    {
        return m_rnd;
    }

    const uint32_t* endRnd() const
    {
        return m_rndEnd;
    }


    uint32_t& scalarState(uint32_t scalarIndex)
    {
        return m_state[scalarIndex * s_regLenWords];
    }

    uint32_t scalarState(uint32_t scalarIndex) const
    {
        return m_state[scalarIndex * s_regLenWords];
    }

    // initializes the first state with a seed
    void __reinit(uint32_t s)
    {
        const uint32_t mask = uint32_t(1812433253UL);
        uint32_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < s_N; i++)
            prev = scalarState(i) = (mask * (prev ^ (prev >> 30)) + i);
    }

protected:
    void reinitPointers()
    {
        m_pst = m_pstEnd;
        m_prnd = (const uint32_t*)(((uint8_t*)m_rnd) + sizeof(m_rnd));
    }

    // extract one of the interleaved state vectors, shift it left by 31 bits and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint32_t* pstate = m_state;
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
        uint32_t* pstate = m_state;
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

    // initializes m_state[s_N] with a seed
    void reinitMainState(uint32_t s)
    {
        __reinit(s);
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    void reinitMainState(const uint32_t* seeds, uint32_t nSeeds)
    {
        __reinit(uint32_t(19650218));
        uint32_t i = 1, j = 0;
        uint32_t k = (s_N > nSeeds ? s_N : nSeeds);
        for (; k; k--) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * uint32_t(1664525)))
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
    }

    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        if ((reinterpret_cast<intptr_t>(m_prnd) % sizeof(m_rnd)) != 0)
            return *m_prnd++;

        if (m_pst != m_pstEnd)
            /* do nothing*/; // most likely case first
        else {
            refill();
        }

        temperRefillBlock<true>(m_pst, m_rnd);
        m_prnd = m_rnd;

        return *m_prnd++;
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        if (m_pst != m_pstEnd)
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
        const uint32_t* pst = m_state;
        for (size_t i = 0; i < s_n32InFullState / s_n32InRndCache; ++i, dst += s_n32InRndCache) {
            temperRefillBlock<false>(pst, dst);
        }
    }

    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        if (size_t nAvailInRnd = std::distance<const uint32_t*>(m_prnd, endRnd()); nAvailInRnd < n) {
            std::copy_n(m_prnd, nAvailInRnd, dst);
            n -= nAvailInRnd;
            dst += nAvailInRnd;
        }
        else {
            std::copy_n(m_prnd, n, dst);
            m_prnd += n;
            return;
        }

        // we are now aligned with Block16
        if (size_t nAvailInState = std::distance(m_pst, end()); nAvailInState < n) {
            if (nAvailInState) {
                n -= nAvailInState;
                do {
                    temperRefillBlock<false>(m_pst, dst);
                    dst += s_n32InRndCache;
                    nAvailInState -= s_n32InRndCache;
                } while (nAvailInState);
            }
        }
        else {
            size_t nFullRndBlocks = n / s_n32InRndCache;
            for (size_t i = 0; i < nFullRndBlocks; ++i) {
                temperRefillBlock<false>(m_pst, dst);
                dst += s_n32InRndCache;
            }
            n = n % s_n32InRndCache;
            if (n) {
                temperRefillBlock<true>(m_pst, m_rnd);
                std::copy_n(m_rnd, n, dst);
                m_prnd = m_rnd + n;
            }
            else {
                m_prnd = endRnd();
            }
            return;
        }

        // we are now aligned with the state vector

        size_t nFullStates = n / s_n32InFullState;
        while (nFullStates--) {
            genrand_uint32_stateBlk(dst);
            dst += s_n32InFullState;
            n -= s_n32InFullState;
        }
        m_pst = end();

        if (n > 0) {
            refill();

            while (n > s_n32InRndCache) {
                temperRefillBlock<false>(m_pst, dst);
                n -= s_n32InRndCache;
                dst += s_n32InRndCache;
            }

            if (n) {
                temperRefillBlock<true>(m_pst, m_rnd);
                std::copy_n(m_rnd, n, dst);
                m_prnd = m_rnd + n;
            }
            else {
                m_prnd = endRnd();
            }
        }
        else {
            m_prnd = endRnd();
        }
    }

public:

    // constructors
    VMT19937Base()
        : m_prnd(nullptr)
        , m_rndEnd(m_rnd + s_n32InRndCache)
        , m_pst(nullptr)
        , m_pstEnd(m_state + s_N * s_n32inReg)
    {}
};

#if (VMT19937_STATIC_CONST==1)
template <size_t RegisterBitLen>
const typename VMT19937Base<RegisterBitLen>::RefillCst VMT19937Base<RegisterBitLen>::s_refillCst;
#endif

} // namespace Details
