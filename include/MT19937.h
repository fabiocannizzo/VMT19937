#pragma once

#include "SIMD.h"
#include "jump_matrix.h"

#include <cstdint>
#include <cstddef>
#include <type_traits>

#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
#  define XVMT_PREFETCH(addr) __builtin_prefetch((addr), 0, 3)
#else
#  define XVMT_PREFETCH(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_T0)
#endif

namespace xvmt {
namespace details {


template <typename WordT, size_t NWords>
struct RndCache
{
    static constexpr bool s_enabled = true;

    RndCache() { setEnd(); }

    void setEnd() { m_cur = m_rnd + NWords; }
    void setBegin() { m_cur = begin(); }
    void setAt(size_t pos) { m_cur = begin() + pos; }

    WordT* begin() { return m_rnd; }
    WordT* current() { return m_cur; }
    const WordT* end() const { return m_rnd + NWords; }
    WordT operator*() const { return *m_cur; }
    WordT& operator*() { return *m_cur; }
    RndCache& operator+=(size_t n) { m_cur += n; return *this; }
    RndCache& operator++() { ++m_cur; return *this; }
    RndCache operator++(int) { RndCache tmp = *this; ++m_cur; return tmp; }
    WordT operator[](size_t i) const { return m_rnd[i]; }

    bool isAtEnd() const { return m_cur == end(); }
    size_t nAvailable() const { return std::distance<const WordT*>(m_cur, end()); }

    alignas(64) WordT m_rnd[NWords];
    WordT* m_cur;
};

template <typename WordT>
struct RndCache<WordT, 0>
{
    void setEnd() {}
    bool isAtEnd() const { return true; }
    size_t nAvailable() const { return 0; }
    WordT* begin() { return nullptr; }
    WordT* current() { return nullptr; }
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
//               true:  only genrand_word_blk() is available.
// Params      - parameter struct selecting 32-bit (MT19937Params) or 64-bit (MT19937_64Params).
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16, typename Params = MT19937Params>
class MT19937Base
{
    static constexpr size_t HwBitLen = IsaTraits<Isa>::HwBitLen;
    static_assert(VRegBitLen == 32 || VRegBitLen == 128 || VRegBitLen == 256 || VRegBitLen == 512,
        "VRegBitLen must be a valid SIMD hardware register width (32, 128, 256, or 512)");
    static_assert(VRegBitLen % Params::s_wordSizeBits == 0,
        "VRegBitLen must be a multiple of the MT word size");
    static_assert(!MonoState || VRegBitLen == HwBitLen,
        "MonoState=true requires VRegBitLen == HwBitLen");

public:
    using word_t = typename Params::word_t;

    static constexpr size_t s_regLenBits = VRegBitLen;
    static constexpr size_t s_regLenBitsHw = HwBitLen;
    static constexpr ISA s_isa = Isa;
    static constexpr int    s_N = Params::s_N;
    static constexpr int    s_M = Params::s_M;
    static constexpr size_t s_nStates = MonoState ? 1 : VRegBitLen / Params::s_wordSizeBits;
    static constexpr size_t s_n32inReg = VRegBitLen / 32;                    // uint32 lanes per logical SIMD register (32-bit path)
    static constexpr size_t s_n32InOneWord = Params::s_n32InOneWord;
    static constexpr size_t s_n32InOneState = Params::s_n32InOneState;
    static constexpr size_t s_n32InFullState = s_n32InOneState * s_nStates;
    static constexpr size_t s_nMatrixBits = Params::s_nMatrixBits;

    using matrix_t = MT19937Matrix;

private:
    static constexpr uint32_t s_cacheLineBytes = 64;
    static_assert(s_cacheLineBytes * 8 >= VRegBitLen, "Assume that the register size is <= than the cache line");
    static constexpr size_t s_n32InBlock = s_cacheLineBytes / sizeof(uint32_t);      // 16 - uint32 slots per cache line
    static constexpr size_t s_nWordsInBlock = s_n32InBlock / Params::s_n32InOneWord; // words per cache line (16 for 32-bit, 8 for 64-bit)
    static_assert(s_n32InOneState * s_nStates % s_n32InBlock == 0, "full state size not divisible by cache line");

    using XV = SimdRegister<s_regLenBits, Isa>;

    [[no_unique_address]] RndCache<word_t, (QryBlk16 ? 0 : s_nWordsInBlock)> m_rndCache;

    const word_t* m_pst;
    const word_t* const m_pstEnd;

protected:
    alignas(64) word_t m_state[s_N * s_nStates];

private:

    // SIMD tempering constants (32-bit path only; cast to uint32_t so they compile for any Params)
    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask1(uint32_t(Params::s_b)), m_mask2(uint32_t(Params::s_c)) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    struct RefillCst
    {
        using XVI = SimdRegister<VRegBitLen, Isa>;
        RefillCst() : m_upperMask(uint32_t(Params::s_upperMask)), m_lowerMask(uint32_t(Params::s_lowerMask)), m_matrixA(uint32_t(Params::s_matrixA)) {}
        const XVI m_upperMask;
        const XVI m_lowerMask;
        const XVI m_matrixA;
    };

    alignas(64) inline static const TemperCst<SimdRegister<s_n32InBlock * 32, Isa>> s_temperCst{};
    alignas(64) inline static const RefillCst s_refillMasks{};

    // SIMD temper (32-bit MT, operates on whole SIMD register)
    template <typename XVI, typename M>
    static FORCE_INLINE XVI temper(XVI y, const M& masks)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & masks.m_mask1);
        y = y ^ ((y << 15) & masks.m_mask2);
        y = y ^ (y >> 18);
        return y;
    }

    // Scalar temper (works for both 32-bit and 64-bit Params)
    static FORCE_INLINE word_t scalarTemper(word_t y)
    {
        y ^= (y >> Params::s_u) & Params::s_d;
        y ^= (y << Params::s_s) & Params::s_b;
        y ^= (y << Params::s_t) & Params::s_c;
        y ^= (y >> Params::s_l);
        return y;
    }

    template <bool Aligned>
    FORCE_INLINE void temperBlock(word_t* dst)
    {
        if constexpr (Params::s_wordSizeBits == 32) {
            static_assert(s_n32InBlock == 16);
            using XVline = SimdRegister<s_n32InBlock * 32, Isa>;
            XVline tmp = temper(XVline(m_pst), s_temperCst);
            tmp.template store<Aligned>(dst);
            m_pst += s_n32InBlock;
        } else {
            for (size_t i = 0; i < s_nWordsInBlock; ++i)
                dst[i] = scalarTemper(m_pst[i]);
            m_pst += s_nWordsInBlock;
        }
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
                    XVMT_PREFETCH(p + n32PerBlk + J1 * s_n32inReg);
                    XVMT_PREFETCH(p + n32PerBlk + JM * s_n32inReg);
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
        if constexpr (Params::s_wordSizeBits == 32) {
            // 32-bit SIMD path (multi-state and mono-state)
            uint32_t* stCur = reinterpret_cast<uint32_t*>(m_state);

            constexpr int N = s_N;
            constexpr int M = s_M;
            static_assert(N == 624 && M == 397, "SIMD unrolling designed for MT19937-32 parameters");

            XV x0(stCur);

            if constexpr (!MonoState) {
                constexpr size_t nUnroll = 4;

                constexpr size_t n1 = (N - M) / nUnroll;
                constexpr size_t r1 = (N - M) % nUnroll;
                stCur = advanceLoop<1, M>(n1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, nUnroll>{});
                if constexpr (r1 > 0)
                    stCur = advanceLoop<1, M>(1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, (int)r1>{});

                constexpr size_t n2 = (M - 1) / nUnroll;
                constexpr size_t r2 = (M - 1) % nUnroll;
                stCur = advanceLoop<1, M - N>(n2, stCur, x0, s_refillMasks, std::make_integer_sequence<int, nUnroll>{});
                if constexpr (r2 > 0)
                    stCur = advanceLoop<1, M - N>(1, stCur, x0, s_refillMasks, std::make_integer_sequence<int, (int)r2>{});

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
        } else {
            // 64-bit scalar path (Phase 1; scalar only, MonoState=true implied for now)
            word_t* st = m_state;
            constexpr int N = s_N;
            constexpr int M = s_M;
            static constexpr word_t mag01[2] = {word_t(0), Params::s_matrixA};
            int i;
            for (i = 0; i < N - M; ++i) {
                word_t x = (st[i] & Params::s_upperMask) | (st[i + 1] & Params::s_lowerMask);
                st[i] = st[i + M] ^ (x >> 1) ^ mag01[x & 1];
            }
            for (; i < N - 1; ++i) {
                word_t x = (st[i] & Params::s_upperMask) | (st[i + 1] & Params::s_lowerMask);
                st[i] = st[i + (M - N)] ^ (x >> 1) ^ mag01[x & 1];
            }
            word_t x = (st[N - 1] & Params::s_upperMask) | (st[0] & Params::s_lowerMask);
            st[N - 1] = st[M - 1] ^ (x >> 1) ^ mag01[x & 1];
        }

        m_pst = m_state;
        m_rndCache.setEnd();
    }

    word_t& scalarState(uint32_t scalarIndex)
    {
        return m_state[scalarIndex * s_nStates];
    }

    word_t scalarState(uint32_t scalarIndex) const
    {
        return m_state[scalarIndex * s_nStates];
    }

    void __reinit(word_t s)
    {
        word_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < (uint32_t)s_N; ++i)
            prev = scalarState(i) = (Params::s_initMul * (prev ^ (prev >> Params::s_initShift)) + i);
        reinitPointers();
    }

protected:
    void reinitPointers()
    {
        m_pst = m_pstEnd;
        m_rndCache.setEnd();
    }

    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        static_assert(Params::s_wordSizeBits == 32, "stateToVector only implemented for 32-bit Params");
        const uint32_t* pstate = reinterpret_cast<const uint32_t*>(m_state);
        pdst[0] = pstate[stateIndex] >> 31;
        for (size_t i = 1; i < (size_t)s_N; ++i) {
            uint32_t word = pstate[i * s_nStates + stateIndex];
            pdst[i - 1] |= word << 1;
            pdst[i] = word >> 31;
        }
    }

    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        static_assert(Params::s_wordSizeBits == 32, "vectorToState only implemented for 32-bit Params");
        uint32_t* pstate = reinterpret_cast<uint32_t*>(m_state);
        const uint32_t* pw = psrc;
        pstate[stateIndex] = 0;
        size_t w;
        for (w = 0; w < (size_t)s_N - 1; ++w) {
            uint32_t word = pw[w];
            pstate[w * s_nStates + stateIndex] |= word << 31;
            pstate[(w + 1) * s_nStates + stateIndex] = word >> 1;
        }
        pstate[w * s_nStates + stateIndex] |= pw[w] << 31;
    }

    void reinitMainState(uint32_t s)
    {
        __reinit(word_t(s));
    }

    // 64-bit scalar seed — only enabled when word_t != uint32_t to avoid redefinition
    template<size_t WB = Params::s_wordSizeBits, std::enable_if_t<WB == 64, int> = 0>
    void reinitMainState(uint64_t s)
    {
        __reinit(s);
    }

    void reinitMainState(const uint32_t* seeds, uint32_t nSeeds)
    {
        static_assert(Params::s_wordSizeBits == 32, "uint32 seed array only supported for 32-bit Params");
        __reinit(word_t(Params::s_arrayInitSeed));
        uint32_t i = 1, j = 0;
        uint32_t k = ((uint32_t)s_N > nSeeds ? (uint32_t)s_N : nSeeds);
        for (; k; --k) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * word_t(Params::s_arrayInitMul1)))
                + seeds[j] + j;
            ++i; ++j;
            if (i >= (uint32_t)s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
        }
        for (k = (uint32_t)s_N - 1; k; --k) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> 30)) * word_t(Params::s_arrayInitMul2)))
                - i;
            ++i;
            if (i >= (uint32_t)s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
        }
        scalarState(0) = Params::s_msb;
    }

    // 64-bit array seed — only enabled when word_t != uint32_t to avoid redefinition
    template<size_t WB = Params::s_wordSizeBits, std::enable_if_t<WB == 64, int> = 0>
    void reinitMainState(const uint64_t* seeds, uint32_t nSeeds)
    {
        __reinit(Params::s_arrayInitSeed);
        uint32_t i = 1, j = 0;
        uint32_t k = ((uint32_t)s_N > nSeeds ? (uint32_t)s_N : nSeeds);
        for (; k; --k) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> Params::s_initShift)) * Params::s_arrayInitMul1))
                + seeds[j] + j;
            ++i; ++j;
            if (i >= (uint32_t)s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
        }
        for (k = (uint32_t)s_N - 1; k; --k) {
            scalarState(i) = (scalarState(i) ^ ((scalarState(i - 1) ^ (scalarState(i - 1) >> Params::s_initShift)) * Params::s_arrayInitMul2))
                - i;
            ++i;
            if (i >= (uint32_t)s_N) { scalarState(0) = scalarState(s_N - 1); i = 1; }
        }
        scalarState(0) = Params::s_msb;
    }

    FORCE_INLINE word_t genrand_word()
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
    FORCE_INLINE void __genrand_word_blk(word_t* dst)
    {
        if (m_pst == m_pstEnd) VM19937_UNLIKELY
            refill();
        temperBlock<false>(dst);
    }

protected:
    // generates a random uint32 on [0, 0xffffffff]
    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(!QryBlk16);
        return (uint32_t)genrand_word();
    }

    // generates a random uint64 on [0, 0xffffffffffffffff]
    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    FORCE_INLINE uint64_t genrand_uint64()
    {
        static_assert(!QryBlk16);
        if constexpr (Params::s_wordSizeBits == 64)
            return genrand_word();
        else
            return (uint64_t(genrand_word()) << 32) | genrand_word();
    }

    template <bool B = QryBlk16, std::enable_if_t<B == QryBlk16, int> = 0>
    void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QryBlk16);
        static_assert(Params::s_wordSizeBits == 32, "blk16 uint32 API only for 32-bit Params");
        __genrand_word_blk(reinterpret_cast<word_t*>(dst));
    }

    // generates s_nWordsInBlock words into dst (cache-line-sized block)
    template <bool B = QryBlk16, std::enable_if_t<B == QryBlk16, int> = 0>
    void genrand_word_blk(word_t* dst)
    {
        static_assert(QryBlk16);
        __genrand_word_blk(dst);
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(!QryBlk16);
        static_assert(Params::s_wordSizeBits == 32, "uint32 anySize API only for 32-bit Params");
        word_t* wdst = reinterpret_cast<word_t*>(dst);
        size_t fromCache = std::min(n, m_rndCache.nAvailable());
        std::copy_n(m_rndCache.current(), fromCache, wdst);
        wdst += fromCache;
        m_rndCache += fromCache;
        n -= fromCache;

        while (n >= s_nWordsInBlock) {
            __genrand_word_blk(wdst);
            wdst += s_nWordsInBlock;
            n -= s_nWordsInBlock;
        }

        if (n > 0) {
            __genrand_word_blk(m_rndCache.begin());
            std::copy_n(m_rndCache.begin(), n, wdst);
            m_rndCache.setAt(n);
        }
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    void genrand_word_anySize(word_t* dst, size_t n)
    {
        static_assert(!QryBlk16);

        size_t fromCache = std::min(n, m_rndCache.nAvailable());
        std::copy_n(m_rndCache.current(), fromCache, dst);
        dst += fromCache;
        m_rndCache += fromCache;
        n -= fromCache;

        while (n >= s_nWordsInBlock) {
            __genrand_word_blk(dst);
            dst += s_nWordsInBlock;
            n -= s_nWordsInBlock;
        }

        if (n > 0) {
            __genrand_word_blk(m_rndCache.begin());
            std::copy_n(m_rndCache.begin(), n, dst);
            m_rndCache.setAt(n);
        }
    }

public:

    MT19937Base()
        : m_pst(nullptr)
        , m_pstEnd(m_state + s_N * s_nStates)
    {
    }
};


} // namespace details
} // namespace xvmt
