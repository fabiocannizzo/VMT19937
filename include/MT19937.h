#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>

namespace xvmt {
namespace details {


template <size_t N32>
struct RndCache
{
    static constexpr bool s_enabled = true;

    RndCache() { setEnd(); }

    void setEnd() { m_cur = m_rnd + N32; }
    void setBegin() { m_cur = begin(); }
    void setAt(size_t pos) { m_cur = begin() + pos; }

    uint32_t *begin() { return m_rnd; }
    uint32_t* current() { return m_cur; }
    const uint32_t* end() const { return m_rnd + N32; }
    uint32_t operator*() const { return *m_cur; }
    uint32_t& operator*() { return *m_cur; }
    RndCache& operator+=(size_t n) { m_cur += n; return *this; }
    RndCache& operator++() { ++m_cur; return *this; }
    RndCache operator++(int) { RndCache tmp = *this; ++m_cur; return tmp; }
    uint32_t operator[](size_t i) const { return m_rnd[i]; }

    bool isAtEnd() const { return m_cur == end(); }
    size_t nAvailable() const { return std::distance<const uint32_t*>(m_cur, end()); }

    alignas(64) uint32_t m_rnd[N32]; // buffer of tempered numbers
    uint32_t* m_cur;
};

template <>
struct RndCache<0>
{
    void setEnd() {}
};

// VRegBitLen  - virtual (logical) SIMD register width in bits. Three constraints apply:
//   (1) Must be a multiple of IsaTraits<Isa>::HwBitLen (enforced by SimdRegister asserts).
//   (2) Must be a multiple of s_wordSizeBits (= 32 for MT).
//   (3) Must correspond to the HwBitLen of a real ISA: one of 32, 128, 256, or 512.
//       This parameter exists for portability: VMT19937<256, ISA::SSE42> and
//       VMT19937<256, ISA::AVX2> produce identical sequences. On SSE42 hardware,
//       each logical 256-bit operation is emulated by two 128-bit hardware instructions.
//       VRegBitLen > HwBitLen is valid; VRegBitLen < HwBitLen is not.
//   For MonoState=true: VRegBitLen must equal HwBitLen (enforced by static_assert below).
// Isa         - target ISA; selects hardware intrinsics and determines HwBitLen.
// MonoState   - false: multi-state vectorized generator (VMT family, nStates > 1);
//               true:  single-state generator using SIMD for intra-state speed (XMT).
// QryBlk16    - false: scalar and any-size query interface enabled;
//               true:  only genrand_uint32_blk16() is available.
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16>
class MT19937Base : public MT19937Params
{
    static constexpr size_t HwBitLen = IsaTraits<Isa>::HwBitLen;
    static_assert(VRegBitLen == 32 || VRegBitLen == 128 || VRegBitLen == 256 || VRegBitLen == 512,
        "VRegBitLen must be a valid SIMD hardware register width (32, 128, 256, or 512)");
    static_assert(VRegBitLen % s_wordSizeBits == 0,
        "VRegBitLen must be a multiple of the MT word size (32)");
    static_assert(!MonoState || VRegBitLen == HwBitLen,
        "MonoState=true requires VRegBitLen == HwBitLen");

public:
    static constexpr size_t s_regLenBits = VRegBitLen;                              // logical SIMD width driving vectorisation (may exceed hardware width)
    static constexpr size_t s_regLenBitsHw = HwBitLen;                         // actual hardware SIMD register width in bits
    static constexpr ISA s_isa = Isa;                                                   // target ISA used for SIMD intrinsic selection
    static constexpr size_t s_nStates = MonoState ? 1 : VRegBitLen / s_wordSizeBits; // parallel MT states packed per logical SIMD register
    static constexpr size_t s_n32inReg = VRegBitLen / 32;                          // uint32 lanes per logical SIMD register
    static constexpr size_t s_n32InFullState = s_n32InOneState * s_nStates;            // 624 * nStates - total uint32 elements in the interleaved state array

    using matrix_t = MT19937Matrix;

private:
    static constexpr uint32_t s_cacheLineBytes = 64;                                    // assumed cache line size in bytes; state array is aligned to this
    static_assert(s_cacheLineBytes * 8 >= VRegBitLen, "Assume that the register size is <= than the cache line");
    static constexpr uint32_t s_n32InBlock = s_cacheLineBytes / sizeof(uint32_t);      // 16 - uint32 elements per cache-line block (one temperBlock call)
    static_assert(s_n32InFullState % s_n32InBlock == 0, "full state size not divisible by cache size");

    using XV = SimdRegister<s_regLenBits, Isa>;

    // This data members is necessary only if QueryMode==QM_Scalar
    [[no_unique_address]] RndCache<(QryBlk16 ? 0 : s_n32InBlock)> m_rndCache; // buffer of tempered numbers

    // This data members are redundant if QueryMode==QM_StateSize
    const uint32_t*m_pst, * const m_pstEnd;    // m_pos==m_pstEnd means the state vector has been consumed and need to be regenerated

protected:
    alignas(64) uint32_t m_state[s_N * s_n32inReg];    // the array of state vectors

private:

    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask1(s_temperMask1), m_mask2(s_temperMask2) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    struct RefillCst
    {
        using XVI = SimdRegister<VRegBitLen, Isa>;
        RefillCst() : m_upperMask(s_upperMask), m_lowerMask(s_lowerMask), m_matrixA(s_matrixA) {}
        const XVI m_upperMask;
        const XVI m_lowerMask;
        const XVI m_matrixA;
    };

    alignas(64) inline static const TemperCst<SimdRegister<s_n32InBlock * 32, Isa>> s_temperCst{};
    alignas(64) inline static const RefillCst s_refillMasks{};

    template <typename XVI, typename M>
    static FORCE_INLINE XVI temper(XVI y, const M& masks)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & masks.m_mask1);
        y = y ^ ((y << 15) & masks.m_mask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <bool Aligned>
    FORCE_INLINE void temperBlock(uint32_t* dst)
    {
        static_assert(s_n32InBlock == 16);

        // a virtual register having the same size as a cache line
        // this may be larger than what available in hardware
        // it is implemented iteratin on the available hardware registers
        using XVline = SimdRegister<s_n32InBlock * 32, Isa>;

        XVline tmp = temper(XVline(m_pst), s_temperCst);
        tmp.template store<Aligned>(dst);

        m_pst += s_n32InBlock;
    }

    static FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm, const RefillCst& masks)
    {
        XV y = XV::bitwiseSelect(masks.m_upperMask, s, sp);
        XV r = sm ^ (y >> 1);
        return r.xorIfOddCst32(sp, masks.m_matrixA);
    }

    template <int J0, int J1, int JM>
    static FORCE_INLINE void multiStateIteration(uint32_t* p, XV& x0, const RefillCst& masks)
    {
        static_assert(!MonoState);
        XV x1(p + J1 * s_n32inReg);
        XV xM(p + JM * s_n32inReg);
        XV tmp = advance1(x0, x1, xM, masks);
        tmp.template store<true>(p + J0 * s_n32inReg);
        x0 = x1;
    }

    template <int J0, int J1, int JM>
    static FORCE_INLINE void monoStateIteration(uint32_t* p, XV& x0, XV& xMlo, const RefillCst& masks)
    {
        static_assert(MonoState);
        constexpr int n32 = (int)s_n32inReg;
        if (n32 == 1) {
            XV xP(p + J1);
            XV xM(p + JM);
            XV r = advance1(x0, xP, xM, masks);
            r.template store<true>(p + J0 * s_n32inReg);
            x0 = xP;
        }
        else {
            constexpr int x1Offset = ((J1 / n32) + (J1 > 0)) * n32;
            XV x1(p + x1Offset);
            constexpr int xMOffset = ((JM / n32) + (JM > 0)) * n32;
            XV xMhi(p + xMOffset);
            XV xP = XV::template alignr32<1>(x0, x1);
            XV xM = XV::template alignr32<s_M % s_n32inReg>(xMlo, xMhi);
            XV r = advance1(x0, xP, xM, masks);
            r.template store<true>(p + J0 * s_n32inReg);
            x0 = x1;
            xMlo = xMhi;
        }
    }

    template <int J1, int JM, int...Is>
    static FORCE_INLINE uint32_t* advanceLoop(size_t nBlkIter, uint32_t* p, XV& x0, const RefillCst& masks, std::integer_sequence<int, Is...>&&)
    {
    	static_assert(!MonoState);
        constexpr size_t nIterPerBlk = sizeof...(Is);
        if constexpr (nIterPerBlk) {
            constexpr size_t n32PerBlk = nIterPerBlk * s_n32inReg;
            auto pend = p + nBlkIter * n32PerBlk;
            do {
                if constexpr (nIterPerBlk >= 4) {
                    __builtin_prefetch(p + n32PerBlk + J1 * s_n32inReg, 0, 3);
                    __builtin_prefetch(p + n32PerBlk + JM * s_n32inReg, 0, 3);
                }
                (multiStateIteration<Is, J1 + Is, JM + Is>(p, x0, masks), ...);
                p += n32PerBlk;
            } while (p != pend);
            return pend;
        }
        else
            return p;
    }

    void NO_INLINE refill()
    {
        uint32_t* stCur = m_state;

        constexpr int N = s_N;
        constexpr int M = s_M;
        static_assert(N == 624 && M == 397, "unrolling designed for these parameters");

        XV x0(stCur);

        if constexpr (!MonoState) {
            constexpr size_t nUnroll = 4;

            // unroll first part of the loop (N-M) iterations
            constexpr size_t n1 = (N - M) / nUnroll;
            constexpr size_t r1 = (N - M) % nUnroll;
            stCur = advanceLoop<1, M>(n1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, nUnroll>{});
            if constexpr (r1 > 0)
                stCur = advanceLoop<1, M>(1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, (int)r1>{});

            // unroll second part of the loop (M-1) iterations
            constexpr size_t n2 = (M - 1) / nUnroll;
            constexpr size_t r2 = (M - 1) % nUnroll;
            stCur = advanceLoop<1, M - N>(n2, stCur, x0, s_refillMasks, std::make_integer_sequence<int, nUnroll>{});
            if constexpr (r2 > 0)
                stCur = advanceLoop<1, M - N>(1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, (int)r2>{});

            // last iteration
            advanceLoop<1 - N, M - N>(1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, 1>{});
        }
        else {
            XV XMlo(stCur + (s_M / s_n32inReg) * s_n32inReg);
            constexpr size_t nIter1 = (N - M) / s_n32inReg;
            for (size_t i = 0; i < nIter1; ++i)
                monoStateIteration<0, 1, M>(stCur + i * s_n32inReg, x0, XMlo, s_refillMasks);
            stCur += nIter1 * s_n32inReg;
            constexpr size_t nIter2 = (M - 1) / s_n32inReg;
            for (size_t i = 0; i < nIter2; ++i)
                monoStateIteration<0, 1, M - N>(stCur + i * s_n32inReg, x0, XMlo, s_refillMasks);
            stCur += nIter2 * s_n32inReg;
            monoStateIteration<0, 1 - N, M - N>(stCur, x0, XMlo, s_refillMasks);
        }

        m_pst = m_state;
        m_rndCache.setEnd();
    }

    uint32_t& scalarState(uint32_t scalarIndex)
    {
        return m_state[scalarIndex * s_nStates];
    }

    uint32_t scalarState(uint32_t scalarIndex) const
    {
        return m_state[scalarIndex * s_nStates];
    }

    // initializes the first state with a seed
    void __reinit(uint32_t s)
    {
        constexpr uint32_t mask = uint32_t(1812433253UL);
        uint32_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < s_N; i++)
            prev = scalarState(i) = (mask * (prev ^ (prev >> 30)) + i);
        reinitPointers();
    }

protected:
    void reinitPointers()
    {
        m_pst = m_pstEnd;
        m_rndCache.setEnd();
    }

    // extract one of the interleaved state vectors, shift it left by 31 bits and save it to dst
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint32_t* pstate = m_state;
        pdst[0] = pstate[stateIndex] >> 31;
        for (size_t i = 1; i < s_N; ++i) {
            uint32_t word = pstate[i * s_nStates + stateIndex];
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
            pstate[w * s_nStates + stateIndex] |= word << 31;
            pstate[(w + 1) * s_nStates + stateIndex] = word >> 1;
        }
        pstate[w * s_nStates + stateIndex] |= pw[w] << 31;
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
    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(!QryBlk16);

        if (!m_rndCache.isAtEnd())
            return *m_rndCache++;

        if (m_pst == m_pstEnd) VM19937_UNLIKELY
            refill();

        temperBlock<true>(m_rndCache.begin());
        m_rndCache.setBegin();

        return *m_rndCache++;
    }

private:
    FORCE_INLINE void __genrand_uint32_blk16(uint32_t* dst)
    {
        if (m_pst == m_pstEnd) VM19937_UNLIKELY
            refill();
        temperBlock<false>(dst);
    }

protected:
    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    template <bool B = QryBlk16, std::enable_if_t<B == QryBlk16, int> = 0>
    void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QryBlk16);
        __genrand_uint32_blk16(dst);
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(!QryBlk16);

        size_t fromCache = std::min(n, m_rndCache.nAvailable());
        std::copy_n(m_rndCache.current(), fromCache, dst);
        dst += fromCache;
        m_rndCache += fromCache;
        n -= fromCache;

        while (n >= 16) {
            __genrand_uint32_blk16(dst);
            dst += 16;
            n -= 16;
        }

        if (n > 0) {
            __genrand_uint32_blk16(m_rndCache.begin());
            std::copy_n(m_rndCache.begin(), n, dst);
            m_rndCache.setAt(n);
        }
    }

public:

    // constructors
    MT19937Base()
        : m_pst(nullptr)
        , m_pstEnd(m_state + s_N * s_nStates)
    {
    }
};


} // namespace details
} // namespace xvmt
