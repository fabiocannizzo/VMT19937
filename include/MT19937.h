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
//   (2) Must be a multiple of s_stateWordBits (= 32 for MT).
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
// Params      - parameter struct: MT19937Params<32> (32-bit) or MT19937Params<64> (64-bit).
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16, typename Params = MT19937Params<32>>
class MT19937Base;


// ============================================================
//  32-bit tempering policy — SIMD block temper
// ============================================================
template <ISA Isa>
struct MT32TemperPolicy
{
private:
    using Params = MT19937Params<32>;

    template <typename XVI, typename M>
    static FORCE_INLINE XVI temper(XVI y, const M& masks)
    {
        y = y ^ (y >> 11);
        y = y ^ ((y << 7) & masks.m_mask1);
        y = y ^ ((y << 15) & masks.m_mask2);
        y = y ^ (y >> 18);
        return y;
    }

    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask1(uint32_t(Params::s_b)), m_mask2(uint32_t(Params::s_c)) {}
        const XVI m_mask1;
        const XVI m_mask2;
    };

    static constexpr size_t s_n32InBlock = 64 / sizeof(uint32_t);  // 16
    using XVline = SimdRegister<s_n32InBlock * 32, Isa>;
    alignas(64) inline static const TemperCst<XVline> s_temperCst{};

public:
    template <bool Aligned, typename Base>
    static FORCE_INLINE void execute(Base& b, typename Base::output_word_t* dst)
    {
        static_assert(s_n32InBlock == 16);
        XVline tmp = temper(XVline(b.m_pst), s_temperCst);
        tmp.template store<Aligned>(dst);
        b.m_pst += s_n32InBlock;
    }
};


// ============================================================
//  32-bit refill policy — SIMD multi/mono-state advance
// ============================================================
template <size_t VRegBitLen, ISA Isa, bool MonoState>
struct MT32RefillPolicy
{
private:
    using Params = MT19937Params<32>;
    using XV = SimdRegister<VRegBitLen, Isa>;
    static constexpr size_t s_n32inReg = VRegBitLen / 32;

    struct RefillCst
    {
        RefillCst() : m_upperMask(uint32_t(Params::s_upperMask)),
                      m_lowerMask(uint32_t(Params::s_lowerMask)),
                      m_matrixA(uint32_t(Params::s_matrixA)) {}
        const XV m_upperMask;
        const XV m_lowerMask;
        const XV m_matrixA;
    };
    alignas(64) inline static const RefillCst s_refillMasks{};

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
            XV xM = XV::template alignr32<Params::s_M % (int)s_n32inReg>(xMlo, xMhi);
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

public:
    template <typename Base>
    static void execute(Base& b)
    {
        uint32_t* stCur = reinterpret_cast<uint32_t*>(b.m_state);

        constexpr int N = Params::s_N;
        constexpr int M = Params::s_M;
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
            XV XMlo(stCur + (M / s_n32inReg) * s_n32inReg);
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

        b.m_pst = b.m_state;
        b.m_rndCache.setEnd();
    }
};


// ============================================================
//  64-bit tempering policy — SIMD block temper
// ============================================================
template <ISA Isa>
struct MT64TemperPolicy
{
private:
    using Params = MT19937Params<64>;

    template <typename XVI>
    static FORCE_INLINE XVI temper(XVI y, const XVI& mask_d, const XVI& mask_b, const XVI& mask_c)
    {
        y = y ^ (shr64(y, Params::s_u) & mask_d);
        y = y ^ (shl64(y, Params::s_s) & mask_b);
        y = y ^ (shl64(y, Params::s_t) & mask_c);
        y = y ^ shr64(y, Params::s_l);
        return y;
    }

    template <typename XVI>
    struct TemperCst
    {
        TemperCst() : m_mask_d(uint64_t(Params::s_d)), m_mask_b(uint64_t(Params::s_b)), m_mask_c(uint64_t(Params::s_c)) {}
        const XVI m_mask_d;
        const XVI m_mask_b;
        const XVI m_mask_c;
    };

    static constexpr size_t s_n64InBlock = 64 / sizeof(uint64_t);  // 8
    using XVline = SimdRegister<s_n64InBlock * 64, Isa>;
    alignas(64) inline static const TemperCst<XVline> s_temperCst{};

public:
    template <bool Aligned, typename Base>
    static FORCE_INLINE void execute(Base& b, typename Base::output_word_t* dst)
    {
        static_assert(s_n64InBlock == 8);
        const auto* pst = reinterpret_cast<const uint32_t*>(b.m_pst);
        XVline tmp = temper(XVline(pst), s_temperCst.m_mask_d, s_temperCst.m_mask_b, s_temperCst.m_mask_c);
        tmp.template store<Aligned>(reinterpret_cast<uint32_t*>(dst));
        b.m_pst += s_n64InBlock;
    }
};


// ============================================================
//  64-bit refill policy — SIMD multi-state / scalar mono-state
// ============================================================
template <size_t VRegBitLen, ISA Isa, bool MonoState>
struct MT64RefillPolicy
{
private:
    using Params = MT19937Params<64>;
    using XV = SimdRegister<VRegBitLen, Isa>;
    static constexpr size_t s_n64inReg = VRegBitLen / 64;
    static constexpr size_t s_n32inReg = VRegBitLen / 32;

    struct RefillCst
    {
        RefillCst() : m_upperMask(uint64_t(Params::s_upperMask)),
                      m_matrixA  (uint64_t(Params::s_matrixA)) {}
        const XV m_upperMask;
        const XV m_matrixA;
    };
    alignas(64) inline static const RefillCst s_refillMasks{};

    static FORCE_INLINE XV advance1(const XV& s, const XV& sp, const XV& sm)
    {
        XV y = XV::bitwiseSelect(s_refillMasks.m_upperMask, s, sp);
        XV r = sm ^ shr64(y, 1);
        return r.xorIfOddCst64(sp, s_refillMasks.m_matrixA);
    }

    static FORCE_INLINE void scalarRefill(typename Params::output_word_t* state, size_t nStates)
    {
        static constexpr typename Params::output_word_t mag01[2] = {0, Params::s_matrixA};
        int i;
        for (i = 0; i < Params::s_N - Params::s_M; ++i)
            for (size_t s = 0; s < nStates; ++s) {
                auto x = (state[i * nStates + s] & Params::s_upperMask)
                       | (state[(i + 1) * nStates + s] & Params::s_lowerMask);
                state[i * nStates + s] = state[(i + Params::s_M) * nStates + s] ^ (x >> 1) ^ mag01[x & 1];
            }
        for (; i < Params::s_N - 1; ++i)
            for (size_t s = 0; s < nStates; ++s) {
                auto x = (state[i * nStates + s] & Params::s_upperMask)
                       | (state[(i + 1) * nStates + s] & Params::s_lowerMask);
                state[i * nStates + s] = state[(i + Params::s_M - Params::s_N) * nStates + s] ^ (x >> 1) ^ mag01[x & 1];
            }
        for (size_t s = 0; s < nStates; ++s) {
            auto x = (state[(Params::s_N - 1) * nStates + s] & Params::s_upperMask)
                   | (state[s] & Params::s_lowerMask);
            state[(Params::s_N - 1) * nStates + s] = state[(Params::s_M - 1) * nStates + s] ^ (x >> 1) ^ mag01[x & 1];
        }
    }

public:
    template <typename Base>
    static void execute(Base& b)
    {
        constexpr int N = Params::s_N;
        constexpr int M = Params::s_M;

        if constexpr (MonoState) {
            scalarRefill(b.m_state, Base::s_nStates);
        }
        else {
            uint32_t* st = reinterpret_cast<uint32_t*>(b.m_state);
            XV x0(st);

            // Phase 1: i = 0..N-M-1
            for (int i = 0; i < N - M; ++i) {
                uint32_t* p = st + i * s_n32inReg;
                XV x1(p + s_n32inReg);
                XV xM(p + M * s_n32inReg);
                XV r = advance1(x0, x1, xM);
                r.template store<true>(p);
                x0 = x1;
            }

            // Phase 2: i = N-M..N-2
            for (int i = N - M; i < N - 1; ++i) {
                uint32_t* p = st + i * s_n32inReg;
                XV x1(p + s_n32inReg);
                XV xM(st + (i + M - N) * s_n32inReg);
                XV r = advance1(x0, x1, xM);
                r.template store<true>(p);
                x0 = x1;
            }

            // Phase 3: i = N-1
            {
                uint32_t* p = st + (N - 1) * s_n32inReg;
                XV x1(st);
                XV xM(st + (M - 1) * s_n32inReg);
                XV r = advance1(x0, x1, xM);
                r.template store<true>(p);
            }
        }

        b.m_pst = b.m_state;
        b.m_rndCache.setEnd();
    }
};


// ============================================================
//  Common base: state storage, cache, and genrand interface
//  Refiller and Temper inject the algorithm-specific hot paths.
// ============================================================
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16, typename Params,
          typename Refiller, typename Temper>
class MT19937BaseImpl
{
    friend Refiller;
    friend Temper;

    static constexpr size_t HwBitLen = IsaTraits<Isa>::HwBitLen;
    static_assert(VRegBitLen == 32 || VRegBitLen == 128 || VRegBitLen == 256 || VRegBitLen == 512,
        "VRegBitLen must be a valid SIMD hardware register width (32, 128, 256, or 512)");
    static_assert(VRegBitLen % Params::s_stateWordBits == 0,
        "VRegBitLen must be a multiple of the MT word size");
    static_assert(!MonoState || VRegBitLen == HwBitLen,
        "MonoState=true requires VRegBitLen == HwBitLen");

public:
    using output_word_t = typename Params::output_word_t;
    using matrix_t = MT19937Matrix<Params::s_stateWordBits>;

    static constexpr size_t s_regLenBits    = VRegBitLen;
    static constexpr size_t s_regLenBitsHw  = HwBitLen;
    static constexpr ISA    s_isa           = Isa;
    static constexpr int    s_N             = Params::s_N;
    static constexpr int    s_M             = Params::s_M;
    static constexpr size_t s_nStates       = MonoState ? 1 : VRegBitLen / Params::s_stateWordBits;
    static constexpr size_t s_n32inReg      = VRegBitLen / 32;
    static constexpr size_t s_n32InOneWord  = Params::s_n32InOneWord;
    static constexpr size_t s_n32InOneState = Params::s_n32InOneState;
    static constexpr size_t s_n32InFullState = s_n32InOneState * s_nStates;
    static constexpr size_t s_nMatrixBits   = Params::s_nMatrixBits;

    // Exposed for policies and RandGen; kept out of private so derived specialisations
    // can also use them in reinitMainState implementations.
    static constexpr size_t s_n32InBlock    = 64 / sizeof(uint32_t);           // 16
    static constexpr size_t s_nWordsInBlock = s_n32InBlock / Params::s_n32InOneWord;

private:
    static_assert(64 * 8 >= VRegBitLen, "Assume that the register size is <= than the cache line");
    static_assert(s_n32InOneState * s_nStates % s_n32InBlock == 0, "full state size not divisible by cache line");

    [[no_unique_address]] RndCache<output_word_t, (QryBlk16 ? 0 : s_nWordsInBlock)> m_rndCache;
    const output_word_t*       m_pst;
    const output_word_t* const m_pstEnd;

protected:
    alignas(64) output_word_t m_state[s_N * s_nStates];

private:
    void NO_INLINE refill()
    {
        Refiller::execute(*this);
    }

    template <bool Aligned>
    FORCE_INLINE void temperBlock(output_word_t* dst)
    {
        Temper::template execute<Aligned>(*this, dst);
    }

protected:
    output_word_t& scalarState(uint32_t scalarIndex)
    {
        return m_state[scalarIndex * s_nStates];
    }

    output_word_t scalarState(uint32_t scalarIndex) const
    {
        return m_state[scalarIndex * s_nStates];
    }

    void __reinit(output_word_t s)
    {
        output_word_t prev = scalarState(0) = s;
        for (uint32_t i = 1; i < (uint32_t)s_N; ++i)
            prev = scalarState(i) = (Params::s_initMul * (prev ^ (prev >> Params::s_initShift)) + i);
        reinitPointers();
    }

    void reinitPointers()
    {
        m_pst = m_pstEnd;
        m_rndCache.setEnd();
    }

    FORCE_INLINE output_word_t genrand_word()
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
    FORCE_INLINE void __genrand_word_blk(output_word_t* dst)
    {
        if (m_pst == m_pstEnd) VM19937_UNLIKELY
            refill();
        temperBlock<false>(dst);
    }

protected:
    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(!QryBlk16);
        return (uint32_t)genrand_word();
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    FORCE_INLINE uint64_t genrand_uint64()
    {
        static_assert(!QryBlk16);
        if constexpr (sizeof(output_word_t) == 4)
            return (uint64_t(genrand_word()) << 32) | genrand_word();
        else
            return genrand_word();
    }

    template <bool B = QryBlk16, std::enable_if_t<B == QryBlk16, int> = 0>
    void genrand_word_blk(output_word_t* dst)
    {
        static_assert(QryBlk16);
        __genrand_word_blk(dst);
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    void genrand_word_anySize(output_word_t* dst, size_t n)
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
    MT19937BaseImpl()
        : m_pst(nullptr)
        , m_pstEnd(m_state + s_N * s_nStates)
    {
    }
};


// ============================================================
//  32-bit specialisation — wires MT32 policies; adds 32-bit
//  state encoding and MT32-only genrand functions.
// ============================================================
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16>
class MT19937Base<VRegBitLen, Isa, MonoState, QryBlk16, MT19937Params<32>>
    : public MT19937BaseImpl<VRegBitLen, Isa, MonoState, QryBlk16, MT19937Params<32>,
                             MT32RefillPolicy<VRegBitLen, Isa, MonoState>,
                             MT32TemperPolicy<Isa>>
{
public:
    using Params  = MT19937Params<32>;
    using base_t  = MT19937BaseImpl<VRegBitLen, Isa, MonoState, QryBlk16, Params,
                                    MT32RefillPolicy<VRegBitLen, Isa, MonoState>,
                                    MT32TemperPolicy<Isa>>;
    using output_word_t = typename Params::output_word_t;

protected:
    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint32_t* pstate = reinterpret_cast<const uint32_t*>(base_t::m_state);
        pdst[0] = pstate[stateIndex] >> 31;
        for (size_t i = 1; i < (size_t)base_t::s_N; ++i) {
            uint32_t word = pstate[i * base_t::s_nStates + stateIndex];
            pdst[i - 1] |= word << 1;
            pdst[i] = word >> 31;
        }
    }

    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        uint32_t* pstate = reinterpret_cast<uint32_t*>(base_t::m_state);
        pstate[stateIndex] = 0;
        size_t w;
        for (w = 0; w < (size_t)base_t::s_N - 1; ++w) {
            uint32_t word = psrc[w];
            pstate[w * base_t::s_nStates + stateIndex] |= word << 31;
            pstate[(w + 1) * base_t::s_nStates + stateIndex] = word >> 1;
        }
        pstate[w * base_t::s_nStates + stateIndex] |= psrc[w] << 31;
    }

    void reinitMainState(uint32_t s)
    {
        base_t::__reinit(output_word_t(s));
    }

    void reinitMainState(const uint32_t* seeds, uint32_t nSeeds)
    {
        base_t::__reinit(output_word_t(Params::s_arrayInitSeed));
        uint32_t i = 1, j = 0;
        uint32_t k = ((uint32_t)base_t::s_N > nSeeds ? (uint32_t)base_t::s_N : nSeeds);
        for (; k; --k) {
            base_t::scalarState(i) = (base_t::scalarState(i) ^ ((base_t::scalarState(i - 1) ^ (base_t::scalarState(i - 1) >> 30)) * output_word_t(Params::s_arrayInitMul1)))
                + seeds[j] + j;
            ++i; ++j;
            if (i >= (uint32_t)base_t::s_N) { base_t::scalarState(0) = base_t::scalarState(base_t::s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
        }
        for (k = (uint32_t)base_t::s_N - 1; k; --k) {
            base_t::scalarState(i) = (base_t::scalarState(i) ^ ((base_t::scalarState(i - 1) ^ (base_t::scalarState(i - 1) >> 30)) * output_word_t(Params::s_arrayInitMul2)))
                - i;
            ++i;
            if (i >= (uint32_t)base_t::s_N) { base_t::scalarState(0) = base_t::scalarState(base_t::s_N - 1); i = 1; }
        }
        base_t::scalarState(0) = Params::s_msb;
    }

    template <bool B = QryBlk16, std::enable_if_t<B == QryBlk16, int> = 0>
    void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QryBlk16);
        base_t::genrand_word_blk(reinterpret_cast<output_word_t*>(dst));
    }

    template <bool B = !QryBlk16, std::enable_if_t<B == !QryBlk16, int> = 0>
    void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(!QryBlk16);
        base_t::genrand_word_anySize(reinterpret_cast<output_word_t*>(dst), n);
    }
};


// ============================================================
//  64-bit specialisation — wires MT64 policies; adds 64-bit
//  state encoding and uint64-seed reinit overloads.
// ============================================================
template <size_t VRegBitLen, ISA Isa, bool MonoState, bool QryBlk16>
class MT19937Base<VRegBitLen, Isa, MonoState, QryBlk16, MT19937Params<64>>
    : public MT19937BaseImpl<VRegBitLen, Isa, MonoState, QryBlk16, MT19937Params<64>,
                             MT64RefillPolicy<VRegBitLen, Isa, MonoState>,
                             MT64TemperPolicy<Isa>>
{
public:
    using Params = MT19937Params<64>;
    using base_t = MT19937BaseImpl<VRegBitLen, Isa, MonoState, QryBlk16, Params,
                                   MT64RefillPolicy<VRegBitLen, Isa, MonoState>,
                                   MT64TemperPolicy<Isa>>;
    using output_word_t = typename Params::output_word_t;

protected:
    // State bit encoding for 64-bit MT (see MT19937_64Matrix::init1_64 for the canonical reference):
    //   Group k (k=0..N-2) occupies bits 64k..64k+63 of the binary vector:
    //     bits  0..32 = x[k] bits 31..63   (33 UPPER_MASK bits of x[k])
    //     bits 33..63 = x[k+1] bits 0..30  (31 LOWER_MASK bits of x[k+1])
    //   Group N-1: only bits 0..32 (UPPER_MASK of x[N-1]).
    //   Total: 64*(N-1)+33 = 19937 bits.

    void stateToVector(size_t stateIndex, uint32_t* pdst) const
    {
        const uint64_t* pstate = reinterpret_cast<const uint64_t*>(base_t::m_state);
        for (size_t k = 0; k < (size_t)base_t::s_N - 1; ++k) {
            uint64_t xk  = pstate[k * base_t::s_nStates + stateIndex];
            uint64_t xk1 = pstate[(k + 1) * base_t::s_nStates + stateIndex];
            pdst[2 * k]     = (uint32_t)(xk >> 31);
            pdst[2 * k + 1] = (uint32_t)((xk >> 63) | ((xk1 & 0x7FFFFFFFULL) << 1));
        }
        uint64_t xLast = pstate[(base_t::s_N - 1) * base_t::s_nStates + stateIndex];
        pdst[2 * (base_t::s_N - 1)]     = (uint32_t)(xLast >> 31);
        pdst[2 * (base_t::s_N - 1) + 1] = (uint32_t)(xLast >> 63);
    }

    void vectorToState(size_t stateIndex, const uint32_t* psrc)
    {
        uint64_t* pstate = reinterpret_cast<uint64_t*>(base_t::m_state);
        // x[0]: LOWER bits (0..30) are not encoded in the state vector; set to zero
        pstate[0 * base_t::s_nStates + stateIndex] =
            ((uint64_t)psrc[0] << 31) | ((uint64_t)(psrc[1] & 1) << 63);
        for (size_t k = 1; k < (size_t)base_t::s_N; ++k) {
            uint64_t upper = ((uint64_t)psrc[2 * k] << 31) | ((uint64_t)(psrc[2 * k + 1] & 1) << 63);
            uint64_t lower = (uint64_t)(psrc[2 * k - 1] >> 1) & 0x7FFFFFFFULL;
            pstate[k * base_t::s_nStates + stateIndex] = upper | lower;
        }
    }

    void reinitMainState(uint32_t s)
    {
        base_t::__reinit(output_word_t(s));
    }

    void reinitMainState(uint64_t s)
    {
        base_t::__reinit(s);
    }

    void reinitMainState(const uint64_t* seeds, uint32_t nSeeds)
    {
        base_t::__reinit(Params::s_arrayInitSeed);
        uint32_t i = 1, j = 0;
        uint32_t k = ((uint32_t)base_t::s_N > nSeeds ? (uint32_t)base_t::s_N : nSeeds);
        for (; k; --k) {
            base_t::scalarState(i) = (base_t::scalarState(i) ^ ((base_t::scalarState(i - 1) ^ (base_t::scalarState(i - 1) >> Params::s_initShift)) * Params::s_arrayInitMul1))
                + seeds[j] + j;
            ++i; ++j;
            if (i >= (uint32_t)base_t::s_N) { base_t::scalarState(0) = base_t::scalarState(base_t::s_N - 1); i = 1; }
            if (j >= nSeeds) j = 0;
        }
        for (k = (uint32_t)base_t::s_N - 1; k; --k) {
            base_t::scalarState(i) = (base_t::scalarState(i) ^ ((base_t::scalarState(i - 1) ^ (base_t::scalarState(i - 1) >> Params::s_initShift)) * Params::s_arrayInitMul2))
                - i;
            ++i;
            if (i >= (uint32_t)base_t::s_N) { base_t::scalarState(0) = base_t::scalarState(base_t::s_N - 1); i = 1; }
        }
        base_t::scalarState(0) = Params::s_msb;
    }
};


} // namespace details
} // namespace xvmt
