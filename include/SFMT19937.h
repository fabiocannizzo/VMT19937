#pragma once

#include "SIMD.h"
#include "jump_matrix.h"
#include "Params.h"

#include <cstdint>
#include <cstddef>

namespace xvmt {
namespace details {

// VRegBitLen  - virtual (logical) SIMD register width in bits. Three constraints apply:
//   (1) Must be a multiple of IsaTraits<Isa>::HwBitLen (enforced by SimdRegister asserts).
//   (2) Must be a multiple of s_stateWordBits (= 128 for SFMT).
//   (3) Must correspond to the HwBitLen of a real ISA: one of 128, 256, or 512.
//       This parameter exists for portability: VSFMT19937<256, ISA::SSE42> and
//       VSFMT19937<256, ISA::AVX2> produce identical sequences. On SSE42 hardware,
//       each logical 256-bit operation is emulated by two 128-bit hardware instructions.
//       VRegBitLen > HwBitLen is valid; VRegBitLen < HwBitLen is not.
// Isa         - target ISA; selects hardware intrinsics and determines HwBitLen.
template <size_t VRegBitLen, ISA Isa>
class SFMT19937Base : public SFMT19937Params
{
    static constexpr size_t s_hwBitLen = IsaTraits<Isa>::s_hwBitLen;
    static_assert(VRegBitLen == 128 || VRegBitLen == 256 || VRegBitLen == 512,
        "VRegBitLen must be a valid SIMD hardware register width (128, 256, or 512)");
    static_assert(VRegBitLen % s_stateWordBits == 0,
        "VRegBitLen must be a multiple of the SFMT word size (128)");

public:
    using output_word_t = uint32_t;

    static constexpr size_t s_regLenBits = VRegBitLen;                              // logical SIMD width driving vectorisation (may exceed hardware width)
    static constexpr size_t s_regLenBitsHw = s_hwBitLen;                         // actual hardware SIMD register width in bits
    static constexpr ISA s_isa = Isa;                                                   // target ISA used for SIMD intrinsic selection
    static constexpr size_t s_nStates = VRegBitLen / s_stateWordBits;               // parallel SFMT states packed per logical SIMD register
    static constexpr size_t s_n32inReg = VRegBitLen / 32;                          // uint32 lanes per logical SIMD register

    static constexpr size_t s_n32InFullState = s_n32InOneState * s_nStates;            // 624 * nStates - total uint32 elements in the interleaved state array

    using matrix_t = SFMT19937Matrix;

private:
    static constexpr size_t s_regLenWords = s_regLenBits / s_stateWordBits;  // FIXME: review this definition

    using XV = SimdRegister<s_regLenBits, Isa>;

protected:
    alignas(64) uint32_t m_state[s_n32InFullState];    // the array of state vectors
    size_t m_step_idx;

private:
    // This data members is necessary only if QueryMode!=QM_StateSize
    const uint32_t* const m_state_end;
    const uint32_t* m_prnd;

    public:
    // Single 128-bit block step update.
    void step()
    {
        uint32_t* st = m_state;
        constexpr int N = s_N;
        alignas(64) static const XV bMask{SFMT19937Params::s_SFMT_MSK1, SFMT19937Params::s_SFMT_MSK2, SFMT19937Params::s_SFMT_MSK3, SFMT19937Params::s_SFMT_MSK4};
        XV xA = XV::template load<true>(st + m_step_idx * s_n32inReg);
        XV xB = XV::template load<true>(st + ((m_step_idx + s_M) % N) * s_n32inReg);
        XV xC = XV::template load<true>(st + ((m_step_idx + N - 2) % N) * s_n32inReg);
        XV xD = XV::template load<true>(st + ((m_step_idx + N - 1) % N) * s_n32inReg);
        XV res = advance1(xA, xB, xC, xD, bMask);
        res.template store<true>(st + m_step_idx * s_n32inReg);
        m_step_idx = (m_step_idx + 1) % N;
        m_prnd = m_state_end; // Force refill cache on next query if needed
    }

    private:

    template <typename XVCst>
    static FORCE_INLINE XV advance1(const XV& xA, const XV& xB, const XV& xC, const XV& xD, const XVCst& bMask)
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
        return XV::bitwiseXorAnd(z ^ x, y, bMask);
    }

    // xA=srcP[JA], xB=xBBase[JA], write to writeBase[WO+JA].
    // Phase 1: xBBase=srcP+s_M, writeBase=dstP, WO=0.
    // Phase 2: xBBase=writeBase=dstReadP, WO=s_N-s_M.
    // A=true: aligned loads/stores (in-place refill on m_state). A=false: unaligned (user buffer).
    template <bool A, int nIter, int JA, int WO, typename XVCst>
    static FORCE_INLINE void unrollD(const uint32_t* srcP, const uint32_t* xBBase,
                                     uint32_t* writeBase, XV& xC, XV& xD, const XVCst& bMask)
    {
        if constexpr (nIter > 0) {
            XV xA = XV::template load<A>(srcP + JA * s_n32inReg);
            XV xB = XV::template load<A>(xBBase + JA * s_n32inReg);
            XV tmp = advance1(xA, xB, xC, xD, bMask);
            xC = xD;
            xD = tmp;
            tmp.template store<A>(writeBase + (WO + JA) * s_n32inReg);
            unrollD<A, nIter - 1, JA + 1, WO>(srcP, xBBase, writeBase, xC, xD, bMask);
        }
    }

    // XBBaseFromSrc=true  (phase 1): xBBase = srcCur + JBstart*..., writeBase = ptr2
    // XBBaseFromSrc=false (phase 2): xBBase = ptr2 = writeBase
    template <bool A, int nUnroll, int nIter, bool XBBaseFromSrc, int JBstart, int WO, typename XVCst>
    static FORCE_INLINE void advanceLoopD(const uint32_t*& srcCur, uint32_t*& ptr2, XV& xC, XV& xD, const XVCst& bMask)
    {
        const size_t nMainIter = nIter / nUnroll;
        if constexpr (nMainIter) {
            auto pend = srcCur + nMainIter * nUnroll * s_n32inReg;
            do {
                const uint32_t* xBBase = XBBaseFromSrc ? srcCur + JBstart * s_n32inReg : ptr2;
                unrollD<A, nUnroll, 0, WO>(srcCur, xBBase, ptr2, xC, xD, bMask);
                srcCur += nUnroll * s_n32inReg;
                ptr2 += nUnroll * s_n32inReg;
            } while (srcCur != pend);
        }
        const size_t nResIter = nIter % nUnroll;
        if constexpr (nResIter) {
            const uint32_t* xBBase = XBBaseFromSrc ? srcCur + JBstart * s_n32inReg : ptr2;
            unrollD<A, nResIter, 0, WO>(srcCur, xBBase, ptr2, xC, xD, bMask);
            srcCur += nResIter * s_n32inReg;
            ptr2 += nResIter * s_n32inReg;
        }
    }

    // A=true: aligned loads/stores (src/dst == m_state). A=false: unaligned (caller buffer).
    // Phase 1: s_N-s_M iters -- xA=src[i], xB=src[i+s_M], write dst[i]
    // Phase 2: s_M iters    -- xA=src[s_N-s_M+i], xB=dst[i], write dst[s_N-s_M+i]
    template <bool A>
    static NO_INLINE void refillImpl(const uint32_t* src, uint32_t* dst)
    {
        const int s_M = SFMT19937Params::s_M;
        alignas(64) static const XV bMask{SFMT19937Params::s_SFMT_MSK1, SFMT19937Params::s_SFMT_MSK2, SFMT19937Params::s_SFMT_MSK3, SFMT19937Params::s_SFMT_MSK4};

        XV xC = XV::template load<A>(src + (s_N - 2) * s_n32inReg);
        XV xD = XV::template load<A>(src + (s_N - 1) * s_n32inReg);

        constexpr size_t nUnroll = (s_regLenBits <= s_regLenBitsHw) ? 8 : 4;

        // phase 1
        const uint32_t* srcCur = src;
        uint32_t* dstCur = dst;
        advanceLoopD<A, nUnroll, s_N - s_M, true,  s_M, 0        >(srcCur, dstCur,      xC, xD, bMask);

        // phase 2
        uint32_t* dstReadCur = dst;
        advanceLoopD<A, nUnroll, s_M,       false, 0,   s_N - s_M>(srcCur, dstReadCur, xC, xD, bMask);
    }

public:
    FORCE_INLINE void refill()
    {
        refillImpl<true>(m_state, m_state);
        m_prnd = begin();
        m_step_idx = 0;
    }

    const uint32_t* begin() const
    {
        return m_state;
    }

    const uint32_t* end() const
    {
        return m_state + s_n32InFullState;
    }

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    //uint32_t& scalarState(uint32_t stateIndex, uint32_t scalarIndex)
    //{
    //    return ((uint32_t*)m_state)[scalarIndex * s_regLenWords];
    //}

    // return the absolute index in the state vector
    // of the 32-bit word `w32RelIndex` belonging to state `stateIndex`
    static size_t w32AbsIndex(size_t w32RelIndex, size_t stateIndex = 0)
    {
        size_t w128RelIndex = w32RelIndex / s_n32InOneWord;
        size_t w128RelOffset = w32RelIndex % s_n32InOneWord;
        return w128RelIndex * s_n32inReg + s_n32InOneWord * stateIndex + w128RelOffset;
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
        uint32_t* psfmt32 = m_state;

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

    void getBlock16(uint32_t* dst)
    {
        std::copy(m_prnd, m_prnd + 16, dst);
        m_prnd += 16;
    }

public:
    // extract one of the interleaved state vectors and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        // For SFMT, we need to extract from interleaved memory but starting at m_step_idx
        for (size_t i = 0; i < (size_t)s_N; ++i) {
            size_t absIdx = (m_step_idx + i) % s_N;
            for (size_t j = 0; j < (size_t)s_n32InOneWord; ++j) {
                pdst[i * s_n32InOneWord + j] = m_state[absIdx * s_n32inReg + stateIndex * s_n32InOneWord + j];
            }
        }
    }

    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        m_step_idx = 0;
        matrixToCube<s_nStates, s_n32InOneWord, uint32_t>(m_state, psrc, s_N, stateIndex);
    }

protected:

    // Initializes the internal state array with a 32-bit integer seed.
    void reinitMainState(uint32_t seed)
    {
        uint32_t* psfmt32 = begin();

        psfmt32[0] = seed;
        for (size_t i = 1; i < s_n32InOneState; i++) {
            psfmt32[w32AbsIndex(i)] = 1812433253UL * (psfmt32[w32AbsIndex(i - 1)]
                ^ (psfmt32[w32AbsIndex(i - 1)] >> 30))
                + i;
        }
        ensure_period();
    }

    // Initializes the internal state array, with an array of 32-bit integers used as the seeds
    void reinitMainState(const uint32_t* init_key, uint32_t key_length)
    {
        uint32_t i, j;
        uint32_t* psfmt32 = m_state;

        const size_t lag = 11;
        const size_t mid = (s_n32InOneState - lag) / 2;

        for (size_t w = 0; w < s_N; ++w)
            for (size_t j = 0; j < s_n32InOneWord; ++j)
                psfmt32[w * s_n32inReg + j] = uint32_t(0x8b8b8b8b);
        size_t count = (key_length + 1 > s_n32InOneState) ? key_length + 1 : s_n32InOneState;
        uint32_t r = func1(psfmt32[w32AbsIndex(0)] ^ psfmt32[w32AbsIndex(mid)] ^ psfmt32[w32AbsIndex(s_n32InOneState - 1)]);
        psfmt32[w32AbsIndex(mid)] += r;
        r += key_length;
        psfmt32[w32AbsIndex(mid + lag)] += r;
        psfmt32[w32AbsIndex(0)] = r;

        count--;
        for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
            r = func1(psfmt32[w32AbsIndex(i)] ^ psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^ psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] += r;
            r += init_key[j] + i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] += r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }
        for (; j < count; j++) {
            r = func1(psfmt32[w32AbsIndex(i)] ^ psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^ psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] += r;
            r += i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] += r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }
        for (j = 0; j < s_n32InOneState; j++) {
            r = func2(psfmt32[w32AbsIndex(i)] + psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] + psfmt32[w32AbsIndex((i + s_n32InOneState - 1) % s_n32InOneState)]);
            psfmt32[w32AbsIndex((i + mid) % s_n32InOneState)] ^= r;
            r -= i;
            psfmt32[w32AbsIndex((i + mid + lag) % s_n32InOneState)] ^= r;
            psfmt32[w32AbsIndex(i)] = r;
            i = (i + 1) % s_n32InOneState;
        }

        ensure_period();
    }

    void reinitPointers()
    {
        m_step_idx = 0;
        m_prnd = end();
    }

    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        if (m_prnd == end()) [[unlikely]]
            refill();
        return *m_prnd++;
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        if (m_prnd == end()) [[unlikely]]
            refill();
        getBlock16(dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        refill();
        std::copy(begin(), end(), dst);
    }

    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        // 1. Consume available from current buffer
        size_t nAvail = std::distance(m_prnd, end());
        if (nAvail > 0) {
            size_t nChunk = std::min(n, nAvail);
            std::copy_n(m_prnd, nChunk, dst);
            m_prnd += nChunk;
            dst += nChunk;
            n -= nChunk;
        }

        if (n == 0) return;

        // 2. Generate blocks directly into dst if possible
        if (n >= s_n32InFullState) {
            // Use the last block in m_state as the source for the first direct refill
            refillImpl<false>(m_state, dst);
            n -= s_n32InFullState;
            dst += s_n32InFullState;

            while (n >= s_n32InFullState) {
                refillImpl<false>(dst - s_n32InFullState, dst);
                n -= s_n32InFullState;
                dst += s_n32InFullState;
            }

            // Sync m_state with the last block generated for future calls
            std::copy_n(dst - s_n32InFullState, s_n32InFullState, m_state);
            m_prnd = end();
        }

        // 3. Final partial block
        if (n > 0) {
            refill();
            std::copy_n(m_prnd, n, dst);
            m_prnd += n;
        }
    }

public:

    // constructors
    SFMT19937Base()
        : m_state_end(m_state + s_n32InFullState)
        , m_step_idx(0)
        , m_prnd(nullptr)
    {}

    SFMT19937Base(const SFMT19937Base& other)
        : m_step_idx(other.m_step_idx)
        , m_state_end(m_state + s_n32InFullState)
        , m_prnd(other.m_prnd == nullptr ? nullptr :
                 m_state + (other.m_prnd - other.m_state))
    {
        std::copy(other.m_state, other.m_state + s_n32InFullState, m_state);
    }

}; // SFMT19937Base

} // namespace details
} // namespace xvmt

