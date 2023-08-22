#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>

namespace Details {

template <size_t RegisterBitLen>
class VSFMT19937Base
{
    static const size_t s_nBits = 19937;
    static const size_t s_wordSizeBits = 128;

    static_assert(RegisterBitLen >= s_wordSizeBits);

public:
    static const size_t s_regLenBits = RegisterBitLen;
    static const size_t s_nStates = RegisterBitLen / s_wordSizeBits;
    static const size_t s_n32InOneWord = s_wordSizeBits / 32;
    typedef BinaryMatrix<156 * 4 * 32> matrix_t;

private:
    const static size_t s_regLenWords = s_regLenBits / s_wordSizeBits;  // FIXME: review this definition
    typedef SimdRegister<s_regLenBits> XV;

    // Period parameters
    static const size_t s_N = s_nBits / s_wordSizeBits + (s_nBits % s_wordSizeBits != 0);
    static_assert(s_N == 156);
    static const size_t s_M = 122;
    static const size_t s_SFMT_N32 = s_N * 4;
    static const size_t s_n32inReg = sizeof(XV) / sizeof(uint32_t);

    static const uint32_t s_rndBlockSize = 64 / sizeof(uint32_t); // exactly one cache line

    alignas(64) XV m_state[s_N];    // the array of state vectors


    // This data members is necessary only if QueryMode!=QM_StateSize
    const uint32_t* m_prnd;

public:
    const static size_t s_qryStateSize = sizeof(m_state) / (sizeof(uint32_t));

private:
    class RefillCst
    {
        static const uint32_t s_SFMT_MSK1 = 0xdfffffefU;
        static const uint32_t s_SFMT_MSK2 = 0xddfecb7fU;
        static const uint32_t s_SFMT_MSK3 = 0xbffaffffU;
        static const uint32_t s_SFMT_MSK4 = 0xbffffff6U;
    public:
        RefillCst() : m_bMask(s_SFMT_MSK1, s_SFMT_MSK2, s_SFMT_MSK3, s_SFMT_MSK4) {}
        const XV m_bMask;
    };


    static FORCE_INLINE XV advance1(const XV& xA, const XV& xB, const XV& xC, const XV& xD, const RefillCst& masks)
    {
        const uint32_t SFMT_SL1 = 18;
        const uint32_t SFMT_SL2 = 1;
        const uint32_t SFMT_SR1 = 11;
        const uint32_t SFMT_SR2 = 1;

        XV y(xB >> SFMT_SR1);
        XV z(XV::template shr128<SFMT_SR2>(xC));
        XV v(xD << SFMT_SL1);
        z = z ^ xA;
        z = z ^ v;
        XV x(XV::template shl128<SFMT_SL2>(xA));
        y = y & masks.m_bMask;
        z = z ^ x;
        return (z ^ y);
    }

    template <int nIter, int JA, int JB>
    static FORCE_INLINE void unroll(XV* p, XV& xC, XV& xD, const RefillCst& masks)
    {
        if constexpr (nIter > 0) {
            XV tmp = advance1(p[JA], p[JB], xC, xD, masks);
            xC = xD;
            xD = tmp;
            p[JA] = tmp;
            unroll<nIter - 1, JA + 1, JB + 1>(p, xC, xD, masks);
        }
    }

    template <int nUnroll, int nIter, int JB>
    static FORCE_INLINE void advanceLoop(XV*& p, XV& xC, XV& xD, const RefillCst& masks)
    {
        if constexpr (nIter >= nUnroll) {
            auto pend = p + (nIter / nUnroll) * nUnroll;
            // unroll the loop in blocks of UnrollBlkSize
            do {
                unroll<nUnroll, 0, JB>(p, xC, xD, masks);
                p += nUnroll;
            } while (p != pend);
        }
        const size_t nResIter = nIter % nUnroll;
        if constexpr (nResIter) {
            unroll<nResIter, 0, JB>(p, xC, xD, masks);
            p += nResIter;
        }
    }

    NO_INLINE void refill()
    {
        static const int N = s_N;
        static const int M = s_M;
        static_assert(N == 156 && M == 122, "unrolling designed for these parameters");

        // Create local copy of the constants and pass them to the function as arguments.
        // Since all functions invoked from here are forced inline, the function arguments
        // will not be passed as arguments via the stack, but reside in CPU registers
        const RefillCst masks{};

        // local variables
        XV* stCur = m_state;
        XV xC = stCur[s_N - 2];
        XV xD = stCur[s_N - 1];

        // unroll first part of the loop: (N-M) iterations
        advanceLoop<2, N - M, M>(stCur, xC, xD, masks);

        // unroll second part of the loop: M iterations
        advanceLoop<2, M, M - N>(stCur, xC, xD, masks);

        m_prnd = begin();
    }

    static size_t w32StateIndex(size_t stateIndex, size_t scalarIndex)
    {
        size_t b32 = scalarIndex / 4;
        size_t o32 = scalarIndex % 4;
        return b32 * s_n32inReg + 4 * stateIndex + o32;
    }

    static size_t idxof(size_t scalarIndex)
    {
        return w32StateIndex(0, scalarIndex);
    }

    const uint32_t* begin() const
    {
        return (const uint32_t*)(m_state);
    }

    const uint32_t* end() const
    {
        return (const uint32_t*)(m_state + s_N);
    }

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    // extract one of the interleaved state vectors and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint32_t* pstate = ((const uint32_t*)m_state) + 4 * stateIndex;
        for (size_t i = 0; i < s_N; ++i)
            for (size_t j = 0; i < 4; ++j)
                pdst[4 * i + j] = pstate[i * s_n32inReg + j];
    }

    // store vector into the interleaved elements of the state vector
    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        uint32_t* pstate = ((const uint32_t*)m_state) + 4 * stateIndex;
        for (size_t i = 0; i < s_N; ++i)
            for (size_t j = 0; j < 4; ++j)
                pstate[i * s_n32inReg + j] = psrc[4 * i + j];
    }


    // This function ensures the period of 2^{MEXP}
    void ensure_period()
    {
        const uint32_t SFMT_PARITY1 = 0x00000001U;
        const uint32_t SFMT_PARITY2 = 0x00000000U;
        const uint32_t SFMT_PARITY3 = 0x00000000U;
        const uint32_t SFMT_PARITY4 = 0x13c9e684U;

        const uint32_t parity[4] = { SFMT_PARITY1, SFMT_PARITY2, SFMT_PARITY3, SFMT_PARITY4 };

        uint32_t inner = 0;
        uint32_t* psfmt32 = (uint32_t*)m_state;

        for (size_t i = 0; i < 4; i++)
            inner ^= psfmt32[i] & parity[i];

        for (size_t i = 16; i > 0; i >>= 1)
            inner ^= inner >> i;

        inner &= 1;

        // check if OK
        if (inner == 1)
            return;

        // check NG, and modification
        for (size_t i = 0; i < 4; i++) {
            uint32_t work = 1;
            for (size_t j = 0; j < 32; j++) {
                if ((work & parity[i]) != 0) {
                    psfmt32[i] ^= work;
                    return;
                }
                work = work << 1;
            }
        }
    }

    // Initializes the internal state array with a 32-bit integer seed.
    void init(uint32_t seed)
    {
        uint32_t* psfmt32 = (uint32_t*)m_state;

        psfmt32[0] = seed;
        for (size_t i = 1; i < s_SFMT_N32; i++) {
            psfmt32[idxof(i)] = 1812433253UL * (psfmt32[idxof(i - 1)]
                ^ (psfmt32[idxof(i - 1)] >> 30))
                + i;
        }
        //sfmt->idx = s_SFMT_N32;
        ensure_period();
    }

    // convenience used in the initialization
    static uint32_t func1(uint32_t x)
    {
        return (x ^ (x >> 27)) * (uint32_t)1664525UL;
    }

    // convenience used in the initialization
    static uint32_t func2(uint32_t x)
    {
        return (x ^ (x >> 27)) * (uint32_t)1566083941UL;
    }

    // Initializes the internal state array, with an array of 32-bit integers used as the seeds
    void init(const uint32_t* init_key, uint32_t key_length)
    {
        uint32_t i, j;
        uint32_t* psfmt32 = (uint32_t*)m_state;

        const size_t lag = 11;
        const size_t mid = (s_SFMT_N32 - lag) / 2;

        memset(m_state, 0x8b, sizeof(m_state));
        size_t count = (key_length + 1 > s_SFMT_N32) ? key_length + 1 : s_SFMT_N32;
        uint32_t r = func1(psfmt32[idxof(0)] ^ psfmt32[idxof(mid)] ^ psfmt32[idxof(s_SFMT_N32 - 1)]);
        psfmt32[idxof(mid)] += r;
        r += key_length;
        psfmt32[idxof(mid + lag)] += r;
        psfmt32[idxof(0)] = r;

        count--;
        for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
            r = func1(psfmt32[idxof(i)] ^ psfmt32[idxof((i + mid) % s_SFMT_N32)] ^ psfmt32[idxof((i + s_SFMT_N32 - 1) % s_SFMT_N32)]);
            psfmt32[idxof((i + mid) % s_SFMT_N32)] += r;
            r += init_key[j] + i;
            psfmt32[idxof((i + mid + lag) % s_SFMT_N32)] += r;
            psfmt32[idxof(i)] = r;
            i = (i + 1) % s_SFMT_N32;
        }
        for (; j < count; j++) {
            r = func1(psfmt32[idxof(i)] ^ psfmt32[idxof((i + mid) % s_SFMT_N32)] ^ psfmt32[idxof((i + s_SFMT_N32 - 1) % s_SFMT_N32)]);
            psfmt32[idxof((i + mid) % s_SFMT_N32)] += r;
            r += i;
            psfmt32[idxof((i + mid + lag) % s_SFMT_N32)] += r;
            psfmt32[idxof(i)] = r;
            i = (i + 1) % s_SFMT_N32;
        }
        for (j = 0; j < s_SFMT_N32; j++) {
            r = func2(psfmt32[idxof(i)] + psfmt32[idxof((i + mid) % s_SFMT_N32)] + psfmt32[idxof((i + s_SFMT_N32 - 1) % s_SFMT_N32)]);
            psfmt32[idxof((i + mid) % s_SFMT_N32)] ^= r;
            r -= i;
            psfmt32[idxof((i + mid + lag) % s_SFMT_N32)] ^= r;
            psfmt32[idxof(i)] = r;
            i = (i + 1) % s_SFMT_N32;
        }

        ensure_period();
    }


    void fillOtherStates(size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        uint32_t* pstate = (uint32_t*)m_state;

#if 0
        // temporary workspace matrix
        BinaryMatrix<2, s_nBits> tmp;

        if (commonJumpRepeat) {
            MYASSERT(commonJump, "commonJump is required when commnJumpRepeat>0");

            // copy state to the first row shifting all bits to the left by 31
            stateToVector(0, (uint32_t*)tmp.rowBegin(0));

            for (size_t i = 0; i < commonJumpRepeat; ++i)
                commonJump->multiplyByColumn(tmp.rowBegin((i + 1) % 2), tmp.rowBegin(i % 2));

            // copy to the state vector shifting all bits to the right by 31
            vectorToState(0, (const uint32_t*)tmp.rowBegin(commonJumpRepeat % 2));
        }
#endif
        if constexpr (s_nStates > 1) {
            if (sequentialJump) {
#if 0
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
#endif
            }
            else {
                for (size_t w = 0; w < s_N; ++w)
                    m_state[w].broadcastLo128();
            }
        }

        m_prnd = end();
    }


    void getBlock16(uint32_t* dst)
    {
        std::copy(m_prnd, m_prnd + 16, dst);
        m_prnd += 16;
    }

protected:
    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        if (m_prnd != end())
            return *m_prnd++;

        refill();

        return *m_prnd++;
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        if (m_prnd != end()) {
            getBlock16(dst);
        }
        else {
            refill();
            getBlock16(dst);
        }
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        refill();
        memcpy(dst, begin(), sizeof(m_state));
    }

public:

    // constructors
    VSFMT19937Base()
        : m_prnd(nullptr)
    {}

    VSFMT19937Base(uint32_t seed, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : VSFMT19937Base()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    VSFMT19937Base(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : VSFMT19937Base()
    {
        reinit(seeds, n_seeds, commonJumpRepeat, commonJump, sequentialJump);
    }

    // initializes m_state[s_N] with a seed
    void reinit(uint32_t s, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        init(s);
        fillOtherStates(commonJumpRepeat, commonJump, sequentialJump);
    }

    // initialize by an array with array-length
    // init_key is the array for initializing keys
    // key_length is its length
    void reinit(const uint32_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        init(seeds, nSeeds);
        fillOtherStates(commonJumpRepeat, commonJump, sequentialJump);
    }

}; // VSFMT19937Base

} // namespace Details

