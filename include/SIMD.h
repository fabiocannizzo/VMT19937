#pragma once

#include "simd_config.h"
#include "portable.h"

#include <algorithm>

namespace xvmt {
namespace details {


// SimdRegister is an abstraction of a packed SIMD register of length VirtualBitLen bits containing words of length 32 bits
// If VirtualBitLen is larger than HwBitLen, then every operation is emulated via iteration over registers of size HwBitLen
// There are specialization for all cases where VirtualBitLen == HwBitLen
template <
    size_t VirtualBitLen,  // length of virtual 32-bit packed register in bits
    ISA Isa                // instruction set architecture
>
struct VirtualRegBase
{
    static constexpr size_t HwBitLen = IsaTraits<Isa>::HwBitLen;
    static_assert(VirtualBitLen >= HwBitLen, "VirtualBitLen must be greater or equal than HwBitLen");
    static_assert(VirtualBitLen % HwBitLen == 0, "VirtualBitLen must be divisble by HwBitLen");
    static_assert(HwBitLen % 32 == 0, "HwBitLen must be divisble by 32");
    static_assert(((HwBitLen/32)& ((HwBitLen/32)-1)) == 0, "HwBitLen/32 must be a power of 2");

    static constexpr size_t s_virtualBitLen = VirtualBitLen;  // virtual register bit length in bits
    static constexpr size_t s_hwBitLen = HwBitLen;            // hardware register bit length in bits
    static constexpr ISA s_isa = Isa;                        // instruction set architecture

protected:
    static_assert(((VirtualBitLen / HwBitLen)& ((VirtualBitLen / HwBitLen)-1)) == 0, "VirtualBitLen / HwBitLen must be a power of 2");
};

// SimdRegister template declaration
template <
    size_t VirtualBitLen,  // length of an abstract virtual register in bits
    ISA Isa,               // instruction set architecture
    typename Enable = void
>
struct SimdRegister;

// SimdRegister: specialization for VirtualBitLen > HwBitLen (emulated via multiple HW registers)
// Excluded: VirtualBitLen==64 is always handled by the explicit SimdRegister<64, Isa, void> below.
template <size_t VirtualBitLen, ISA Isa>
struct SimdRegister<VirtualBitLen, Isa, std::enable_if_t<(VirtualBitLen > IsaTraits<Isa>::HwBitLen && VirtualBitLen != 64), void>>
    : VirtualRegBase<VirtualBitLen, Isa>
{
private:
    static constexpr size_t HwBitLen = VirtualRegBase<VirtualBitLen, Isa>::HwBitLen;
    static constexpr size_t s_M = VirtualBitLen / HwBitLen;
    static constexpr size_t N32 = VirtualBitLen / 32;
    static constexpr size_t N128 = VirtualBitLen / 128;
    using  XVHw = SimdRegister<HwBitLen, Isa>;

    struct Aux
    {
        Aux() : ar{ {} } {}
        Aux(XVHw v) { std::fill_n(ar, s_M, v); }
        Aux(uint32_t v) : Aux(XVHw(v)) {}
        XVHw& operator[](size_t i) { return ar[i]; }
        const XVHw& operator[](size_t i) const { return ar[i]; }
        const XVHw* begin() const { return ar; }
        XVHw* begin() { return ar; }
        const XVHw* end() const { return ar + s_M; }
        XVHw* end() { return ar + s_M; }
    private:
        XVHw ar[s_M];
    };

public:
    Aux m_v;

    using XV = SimdRegister<VirtualBitLen, Isa>;

    SimdRegister() {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(v) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(XVHw(v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        if constexpr (HwBitLen == 32) {
            for (size_t i = 0; i < N128; ++i) {
                m_v[0 + 4 * i] = v0;
                m_v[1 + 4 * i] = v1;
                m_v[2 + 4 * i] = v2;
                m_v[3 + 4 * i] = v3;
            }
        }
        else {
            m_v = Aux(XVHw(v0, v1, v2, v3));
        }
    }
    FORCE_INLINE SimdRegister(const uint32_t* p)
    {
        for (size_t i = 0; i < s_M; ++i, p += sizeof(XVHw) / sizeof(uint32_t))
            m_v[i] = XVHw(p);
    }
    FORCE_INLINE SimdRegister(const Aux& v) : m_v(v) {}

    template <bool A = false>
    FORCE_INLINE void store(uint32_t* p) const
    {
        for (size_t i = 0; i < s_M; ++i, p += sizeof(XVHw) / sizeof(uint32_t))
            m_v[i].template store<A>(p);
    }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p_void)
    {
        const uint32_t* p = (const uint32_t*)p_void;
        XV r;
        for (size_t i = 0; i < s_M; ++i, p += sizeof(XVHw) / sizeof(uint32_t))
            r.m_v[i] = XVHw::template load<A>(p);
        return r;
    }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = a.m_v[i] & b.m_v[i]; return r; }
    template <typename XVI>
    friend FORCE_INLINE XV operator&(const XV& a, const XVI& b)
    {
        XV r;
        if constexpr (XVI::s_virtualBitLen > XVI::s_hwBitLen) {
            const size_t N = XVI::s_virtualBitLen / XVI::s_hwBitLen;
            for (size_t i = 0; i < s_M / N; ++i)
                for (size_t j = 0; j < N; ++j)
                    r.m_v[i * N + j] = a.m_v[i * N + j] & b.m_v[j];
        }
        else {
            for (size_t i = 0; i < s_M; ++i)
                r.m_v[i] = a.m_v[i] & b;
        }
        return r;
    }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = a.m_v[i] ^ b.m_v[i]; return r; }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = a.m_v[i] | b.m_v[i]; return r; }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = a.m_v[i] << n; return r; }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = a.m_v[i] >> n; return r; }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = shl64(a.m_v[i], n); return r; }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { XV r; for (size_t i = 0; i < s_M; ++i) r.m_v[i] = shr64(a.m_v[i], n); return r; }

    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = m_v[i].xorIfOddCst64(cond.m_v[i], cst.m_v[i]);
        return r;
    }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b).
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = XVHw::bitwiseSelect(mask.m_v[i], a.m_v[i], b.m_v[i]);
        return r;
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = XVHw::bitwiseXorAnd(a.m_v[i], b.m_v[i], c.m_v[i]);
        return r;
    }

    // Carry-less multiplication of two 64-bit polynomials.
    // imm8: 0x00: lo-lo, 0x01: hi-lo, 0x10: lo-hi, 0x11: hi-hi
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
        XV r;
        r.m_v[0] = XVHw::template clmul<imm8>(a.m_v[0], b.m_v[0]);
        for (size_t i = 1; i < s_M; ++i)
            r.m_v[i] = XVHw::zero();
        return r;
    }

    // Per-128-bit lane shift left by nBytes.
    template <int nBytes>
    FORCE_INLINE static XV shl128(const XV& a)
    {
        static_assert(N128 > 0);
        XV r;
        if constexpr (HwBitLen == 32) {
            for (size_t s = 0; s < N128; ++s) {
                XVHw current = a.m_v[4 * s];
                r.m_v[4 * s] = current << 8;
                for (size_t j = 1; j < 4; ++j) {
                    XVHw prev = current;
                    current = a.m_v[4 * s + j];
                    r.m_v[4 * s + j] = (current << 8) | (prev >> 24);
                }
            }
        }
        else {
            for (size_t i = 0; i < s_M; ++i)
                r.m_v[i] = XVHw::template shl128<nBytes>(a.m_v[i]);
        }
        return r;
    }

    template <int nBytes>
    FORCE_INLINE static XV shr128(const XV& a)
    {
        static_assert(N128 > 0);
        XV r;
        if constexpr (HwBitLen == 32) {
            for (size_t s = 0; s < N128; ++s) {
                XVHw current = a.m_v[4 * s + 3];
                r.m_v[4 * s + 3] = current >> 8;
                for (int j = 2; j >= 0; --j) {
                    XVHw next = current;
                    current = a.m_v[4 * s + j];
                    r.m_v[4 * s + j] = (current >> 8) | (next << 24);
                }
            }
        }
        else {
            for (size_t i = 0; i < s_M; ++i)
                r.m_v[i] = XVHw::template shr128<nBytes>(a.m_v[i]);
        }
        return r;
    }

    FORCE_INLINE bool eq(const XV& rhs) const
    {
        for (size_t i = 0; i < s_M; ++i)
            if (!m_v[i].eq(rhs.m_v[i])) return false;
        return true;
    }

    static FORCE_INLINE XV zero() { return XV(uint32_t(0)); }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = m_v[i].ifOddCst32ElseZero(cst32.m_v[i]);
        return r;
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = m_v[i].xorIfOddCst32(cond.m_v[i], cst.m_v[i]);
        return r;
    }

    uint8_t parity() const
    {
        uint8_t p = 0;
        for (size_t i = 0; i < s_M; ++i)
            p ^= m_v[i].parity();
        return p;
    }
};

template <ISA Isa>
struct SimdRegister<32, Isa, void>
{
    static const size_t s_virtualBitLen = 32;
    static const size_t s_hwBitLen = 32;
    static const ISA s_isa = Isa;

    uint32_t m_v;

    typedef SimdRegister<32, Isa> XV;

    SimdRegister() : m_v(0) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(*(const uint32_t*)p) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { *dst = m_v; }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { return *(const uint32_t*)p; }

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }

    friend FORCE_INLINE XV operator<<(const XV a, const int n) { return a.m_v << n; }
    friend FORCE_INLINE XV operator>>(const XV a, const int n) { return a.m_v >> n; }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    static FORCE_INLINE XV zero() { return uint32_t(0); }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return (m_v & 1) ? cst.m_v : 0; }
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const { return ifOddCstThenZero(cst32); }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b) { return (mask.m_v & a.m_v) | (~mask.m_v & b.m_v); }
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c) { return a.m_v ^ (b.m_v & c.m_v); }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const { return *this ^ XV(cond.m_v & 1 ? cst.m_v : 0); }

    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 1);
        if constexpr (n32FromSecond == 0) return a;
        else return b;
    }

    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b) {
        uint64_t res = 0;
        for (int i = 0; i < 32; i++) if ((a.m_v >> i) & 1) res ^= (uint64_t(b.m_v) << i);
        return uint32_t(res);
    }

    uint8_t parity() const { return popcnt(m_v) % 2; }
};


template <ISA Isa>
struct SimdRegister<64, Isa, void>
{
    static const size_t s_virtualBitLen = 64;
    static const size_t s_hwBitLen = 64;
    static const ISA s_isa = Isa;

    uint64_t m_v;

    typedef SimdRegister<64, Isa> XV;

    SimdRegister() : m_v(0) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(*(const uint64_t*)p) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { *(uint64_t*)dst = m_v; }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { return *(const uint64_t*)p; }

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }

    friend FORCE_INLINE XV shr64(const XV a, int n) { return uint64_t(a.m_v >> n); }
    friend FORCE_INLINE XV shl64(const XV a, int n) { return uint64_t(a.m_v << n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    static FORCE_INLINE XV zero() { return uint64_t(0); }

    // Conditional masking: result = (this & 1) ? cst : 0.
    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const
    {
        if constexpr (Isa == ISA::AVX2 || Isa == ISA::AVX512) {
            return (m_v & 1) * cst.m_v;
        }
        else {
            return (m_v & 1) ? cst.m_v : uint64_t(0);
        }
    }

    FORCE_INLINE XV ifOddCst64ElseZero(const XV cst64) const { return ifOddCstThenZero(cst64); }

    // Bitwise Selection: result = (mask & a) | (~mask & b).
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return (mask.m_v & a.m_v) | (~mask.m_v & b.m_v);
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
        return a.m_v ^ (b.m_v & c.m_v);
    }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        return *this ^ XV(cond.m_v & 1 ? cst.m_v : uint64_t(0));
    }

    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 2);
        if constexpr (n32FromSecond == 0) return a;
        else if constexpr (n32FromSecond == 2) return b;
        else return (a.m_v >> 32) | (b.m_v << 32);
    }

    // Carry-less multiplication of two 64-bit polynomials (scalar fallback).
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
        uint64_t res_lo = 0;
        for (int i = 0; i < 64; i++) {
            if ((a.m_v >> i) & 1) res_lo ^= (b.m_v << i);
        }
        return res_lo;
    }
};

#if SIMD_N_BITS>=128
#if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__) || defined(_M_ARM64) || defined(__arm__)
template <>
struct SimdRegister<128, ISA::NEON, void> : VirtualRegBase<128, ISA::NEON>
{
    uint32x4_t m_v;

    typedef SimdRegister<128, ISA::NEON> XV;

    SimdRegister() {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(vdupq_n_u32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(vreinterpretq_u32_u64(vdupq_n_u64(v))) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        alignas(16) uint32_t data[4] = { v0, v1, v2, v3 };
        m_v = vld1q_u32(data);
    }
    FORCE_INLINE SimdRegister(const void* p) : m_v(vld1q_u32((const uint32_t*)p)) {}
    FORCE_INLINE SimdRegister(uint32x4_t v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { vst1q_u32(dst, m_v); }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { return vld1q_u32((const uint32_t*)p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return vandq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return veorq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return vorrq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator*(const XV& a, const XV& b) { return vmulq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return vshlq_u32(a.m_v, vdupq_n_s32(n)); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return vshrq_n_u32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return vreinterpretq_u32_u64(vshlq_u64(vreinterpretq_u64_u32(a.m_v), vdupq_n_s64( n))); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return vreinterpretq_u32_u64(vshlq_u64(vreinterpretq_u64_u32(a.m_v), vdupq_n_s64(-n))); }

    FORCE_INLINE bool eq(const XV& rhs) const
    {
        uint32x4_t cmp = vceqq_u32(m_v, rhs.m_v);
        // All bits must be 1 for each lane
        return vgetq_lane_u32(cmp, 0) == 0xFFFFFFFFU &&
               vgetq_lane_u32(cmp, 1) == 0xFFFFFFFFU &&
               vgetq_lane_u32(cmp, 2) == 0xFFFFFFFFU &&
               vgetq_lane_u32(cmp, 3) == 0xFFFFFFFFU;
    }

    // Concatenates registers a and b, then extracts a register-sized window starting from the n32-th word.
    // Effectively shifts the combined [a, b] window left by n32 words.
    // Example (n32FromSecond=1, 128-bit):
    //   a = {a0, a1, a2, a3}
    //   b = {b0, b1, b2, b3}
    //   result = {a1, a2, a3, b0}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 4, "n32FromSecond must be <=4 with NEON");
        if constexpr (n32FromSecond == 0) return a;
        else if constexpr (n32FromSecond == 4) return b;
        else return vextq_u32(a.m_v, b.m_v, n32FromSecond);
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a)
    {
        if constexpr (n == 0) return a;
        else if constexpr (n >= 16) return zero();
        else return vreinterpretq_u32_u8(vextq_u8(vdupq_n_u8(0), vreinterpretq_u8_u32(a.m_v), 16 - n));
    }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a)
    {
        if constexpr (n == 0) return a;
        else if constexpr (n >= 16) return zero();
        else return vreinterpretq_u32_u8(vextq_u8(vreinterpretq_u8_u32(a.m_v), vdupq_n_u8(0), n));
    }

    static FORCE_INLINE XV zero() { return vdupq_n_u32(0); }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b). Uses vbsl instruction.
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return vbslq_u32(mask.m_v, a.m_v, b.m_v);
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
        return a ^ (b & c);
    }

    // Conditional masking: For each 32-bit lane i, result[i] = (this[i] & 1) ? cst32[i] : 0.
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        int32x4_t mask = vshrq_n_s32(vreinterpretq_s32_u32(vshlq_n_u32(m_v, 31)), 31);
        return vandq_u32(vreinterpretq_u32_s32(mask), cst32.m_v);
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    // Conditional masking: For each 64-bit lane i, result[i] = (this[i] & 1) ? cst64[i] : 0.
    FORCE_INLINE XV ifOddCst64ElseZero(const XV cst64) const
    {
        int64x2_t mask = vshrq_n_s64(vreinterpretq_s64_u64(vshlq_n_u64(vreinterpretq_u64_u32(m_v), 63)), 63);
        return vandq_u32(vreinterpretq_u32_s64(mask), cst64.m_v);
    }

    // Conditional XOR: For each 32-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst64ElseZero(cst);
    }

    // Carry-less multiplication of two 64-bit polynomials.
    // imm8: 0x00: lo-lo, 0x01: hi-lo, 0x10: lo-hi, 0x11: hi-hi
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
#if (defined(__ARM_FEATURE_CRYPTO) || defined(__ARM_FEATURE_AES)) && !defined(_MSC_VER)
        uint64_t a64, b64;
        if constexpr ((imm8 & 0x01) == 0) a64 = vgetq_lane_u64(vreinterpretq_u64_u32(a.m_v), 0);
        else a64 = vgetq_lane_u64(vreinterpretq_u64_u32(a.m_v), 1);
        if constexpr ((imm8 & 0x10) == 0) b64 = vgetq_lane_u64(vreinterpretq_u64_u32(b.m_v), 0);
        else b64 = vgetq_lane_u64(vreinterpretq_u64_u32(b.m_v), 1);
        return vreinterpretq_u32_p128(vmull_p64((poly64_t)a64, (poly64_t)b64));
#else
        // Fallback for non-crypto ARM or MSVC
        uint64_t a64 = (imm8 & 0x01) == 0 ? vgetq_lane_u64(vreinterpretq_u64_u32(a.m_v), 0) : vgetq_lane_u64(vreinterpretq_u64_u32(a.m_v), 1);
        uint64_t b64 = (imm8 & 0x10) == 0 ? vgetq_lane_u64(vreinterpretq_u64_u32(b.m_v), 0) : vgetq_lane_u64(vreinterpretq_u64_u32(b.m_v), 1);

        uint64_t res_lo = 0, res_hi = 0;
        for (int i = 0; i < 64; i++) {
            if ((a64 >> i) & 1) {
                res_lo ^= (b64 << i);
                if (i > 0) res_hi ^= (b64 >> (64 - i));
            }
        }
        alignas(16) uint64_t res[2] = { res_lo, res_hi };
        return XV(res);
#endif
    }

    // Returns the parity (reduction XOR) of all bits in the register.
    uint8_t parity() const
    {
        uint32_t d = vgetq_lane_u32(m_v, 0) ^ vgetq_lane_u32(m_v, 1) ^
                     vgetq_lane_u32(m_v, 2) ^ vgetq_lane_u32(m_v, 3);
        return popcnt(d) & 1;
    }
};
#else
template <ISA Isa>
struct SimdRegister<128, Isa, std::enable_if_t<Isa == ISA::SSE2 || Isa == ISA::SSE42 || Isa == ISA::AVX2 || Isa == ISA::AVX512, void>>
{
    static constexpr size_t s_virtualBitLen = 128;
    static constexpr size_t s_hwBitLen = 128; // Logically 128 for this specialization
    static constexpr ISA s_isa = Isa;

    __m128i m_v;

    typedef SimdRegister<128, Isa> XV;

    SimdRegister() : m_v(_mm_undefined_si128()) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(_mm_set1_epi64x((long long)v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm_setr_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm_load_si128((const __m128i*)p)) {}
    FORCE_INLINE SimdRegister(__m128i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { if (A) _mm_store_si128((__m128i*)dst, m_v); else _mm_storeu_si128((__m128i*)dst, m_v); }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { if constexpr (A) return _mm_load_si128((const __m128i*)p); else return _mm_loadu_si128((const __m128i*)p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator*(const XV& a, const XV& b) { return _mm_mullo_epi32(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm_srli_epi64(a.m_v, n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return _mm_test_all_ones(_mm_cmpeq_epi32(m_v, rhs.m_v)); }

    // Concatenates registers a and b, then extracts a register-sized window starting from the n32-th word.
    // Effectively shifts the combined [a, b] window left by n32 words.
    // Example (n32FromSecond=1, 128-bit):
    //   a = {a0, a1, a2, a3}
    //   b = {b0, b1, b2, b3}
    //   result = {a1, a2, a3, b0}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 4, "n32FromSecond must be <=4 with SSE");
        if constexpr (n32FromSecond == 0)
            return a;
        else if constexpr (n32FromSecond == 4)
            return b;
        else
            return _mm_alignr_epi8(b.m_v, a.m_v, 4*n32FromSecond);
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm_bslli_si128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm_bsrli_si128(a.m_v, n); }

    static FORCE_INLINE XV zero() { return _mm_setzero_si128(); }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b).
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
    #if defined(__AVX512VL__)
        return _mm_ternarylogic_epi32(mask.m_v, a.m_v, b.m_v, 0xCA);
    #else
        return b ^ (mask & (a ^ b));
    #endif
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
    #if defined(__AVX512VL__)
        return _mm_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, 0x78);
    #else
        return a ^ (b & c);
    #endif
    }


    // Conditional masking: For each 32-bit lane i, result[i] = (this[i] & 1) ? cst32[i] : 0.
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        const __m128i z = _mm_setzero_si128();
        const __m128i lowestBit = _mm_slli_epi32(m_v, 31);
        const __m128i isOdd = _mm_cmpgt_epi32(z, lowestBit);
        return _mm_and_si128(isOdd, cst32.m_v);
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    // Conditional masking: For each 64-bit lane i, result[i] = (this[i] & 1) ? cst64[i] : 0.
    FORCE_INLINE XV ifOddCst64ElseZero(const XV cst64) const
    {
        // sign-extend bit 0 of each 64-bit lane to a full-lane mask
        __m128i lowestBit = _mm_slli_epi64(m_v, 63);
        __m128i hi   = _mm_shuffle_epi32(lowestBit, 0xF5);  // broadcast high-32 of each 64-bit lane
        __m128i mask = _mm_srai_epi32(hi, 31);
        return _mm_and_si128(mask, cst64.m_v);
    }

    // Conditional XOR: For each 32-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst64ElseZero(cst);
    }

    // Carry-less multiplication of two 64-bit polynomials (PCLMULQDQ).
    // imm8: 0x00: lo-lo, 0x01: hi-lo, 0x10: lo-hi, 0x11: hi-hi
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
        return _mm_clmulepi64_si128(a.m_v, b.m_v, imm8);
    }

    // Returns the parity (reduction XOR) of all bits in the register.
    uint8_t parity() const
    {
        __m128i hi(_mm_shuffle_epi32(m_v, 2 | (3 << 2)));
        __m128i mix = _mm_xor_si128(hi, m_v);
        uint64_t d = _mm_extract_epi64(mix, 0);
        return popcnt(d) & 1;
    }
};
#endif
#endif

#if SIMD_N_BITS>=256

#if defined(__ARM_FEATURE_SVE) && (defined(__aarch64__) || defined(_M_ARM64))
// Fallback for SVE/NEON bridge intrinsics if missing
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ >= 10)
    // Some GCC versions might be missing these or have them under different names
    // We use a local helper to avoid name conflicts.
    static FORCE_INLINE uint32x4_t safe_svget_neonq_u32(svuint32_t v) {
        alignas(16) uint32_t buf[4];
        svbool_t pg4 = svwhilelt_b32(0, 4);
        svst1_u32(pg4, buf, v);
        return vld1q_u32(buf);
    }
    static FORCE_INLINE svuint32_t safe_svset_neonq_u32(svuint32_t v, uint32x4_t neon) {
        alignas(16) uint32_t buf[4];
        vst1q_u32(buf, neon);
        svbool_t pg4 = svwhilelt_b32(0, 4);
        // We want to replace the bottom 128 bits. svld1_u32(pg, p) loads into a new vector with 0s elsewhere.
        // We need to merge it.
        return svsel_u32(pg4, svld1_u32(pg4, buf), v);
    }
    #define svget_neonq_u32 safe_svget_neonq_u32
    #define svset_neonq_u32 safe_svset_neonq_u32
#else
    #define safe_svget_neonq_u32 svget_neonq_u32
    #define safe_svset_neonq_u32 svset_neonq_u32
#endif

#if defined(__GNUC__) && !defined(__clang__) && defined(__ARM_FEATURE_SVE_BITS) && (__ARM_FEATURE_SVE_BITS == 256)
    typedef svuint32_t svuint32_fixed_t __attribute__((arm_sve_vector_bits(256)));
#else
    typedef svuint32_t svuint32_fixed_t;
#endif

template <>
struct SimdRegister<256, ISA::SVE256, void> : VirtualRegBase<256, ISA::SVE256>
{
    svuint32_fixed_t m_v;

    typedef SimdRegister<256, ISA::SVE256> XV;

    SimdRegister() {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(svdup_u32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(svreinterpret_u32_u64(svdup_u64(v))) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        alignas(32) uint32_t data[8] = { v0, v1, v2, v3, v0, v1, v2, v3 };
        m_v = svld1_u32(svptrue_b32(), data);
    }
    FORCE_INLINE SimdRegister(const void* p)
        : m_v(svld1_u32(svptrue_b32(), (const uint32_t*)p)) {}
    FORCE_INLINE SimdRegister(svuint32_t v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { svst1_u32(svptrue_b32(), dst, m_v); }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p)
    {
        return svld1_u32(svptrue_b32(), (const uint32_t*)p);
    }

    static FORCE_INLINE XV zero() { return svdup_u32(0); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b)
    { return svand_u32_x(svptrue_b32(), a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b)
    { return sveor_u32_x(svptrue_b32(), a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b)
    { return svorr_u32_x(svptrue_b32(), a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n)
    { return svlsl_u32_x(svptrue_b32(), a.m_v, svdup_u32(n)); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n)
    { return svlsr_u32_x(svptrue_b32(), a.m_v, svdup_u32(n)); }

    friend FORCE_INLINE XV shl64(const XV& a, int n)
    {
        return svreinterpret_u32_u64(
            svlsl_u64_x(svptrue_b64(), svreinterpret_u64_u32(a.m_v), svdup_u64(n)));
    }
    friend FORCE_INLINE XV shr64(const XV& a, int n)
    {
        return svreinterpret_u32_u64(
            svlsr_u64_x(svptrue_b64(), svreinterpret_u64_u32(a.m_v), svdup_u64(n)));
    }

    // Concatenates registers a and b, then extracts a register-sized window starting from the n32-th word.
    // Effectively shifts the combined [a, b] window left by n32 words.
    // Example (n32FromSecond=1, 256-bit):
    //   a = {a0, a1, a2, a3, a4, a5, a6, a7}
    //   b = {b0, b1, b2, b3, b4, b5, b6, b7}
    //   result = {a1, a2, a3, a4, a5, a6, a7, b0}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 8, "n32FromSecond must be <= 8 for SVE256");
        if constexpr (n32FromSecond == 0) return a;
        if constexpr (n32FromSecond == 8) return b;
        return svext_u32(a.m_v, b.m_v, n32FromSecond);
    }

    // Per-128-bit-lane byte shift -- matches AVX2 _mm256_bslli_epi128 / _mm256_bsrli_epi128 semantics.
    // Each 128-bit lane is shifted independently, using NEON intrinsics on the extracted q-registers.
    template <int n>
    static FORCE_INLINE XV shl128(const XV& a)
    {
        if constexpr (n == 0) return a;
        else if constexpr (n >= 16) return zero();
        svuint32_t z = svdup_u32(0);
        uint32x4_t lo = svget_neonq_u32(a.m_v);
        uint32x4_t hi = svget_neonq_u32(svext_u32(a.m_v, z, 4));
        lo = vreinterpretq_u32_u8(vextq_u8(vdupq_n_u8(0), vreinterpretq_u8_u32(lo), 16 - n));
        hi = vreinterpretq_u32_u8(vextq_u8(vdupq_n_u8(0), vreinterpretq_u8_u32(hi), 16 - n));
        svuint32_t lo_sve = svset_neonq_u32(z, lo);
        svuint32_t hi_sve = svset_neonq_u32(z, hi);
        return svorr_u32_x(svptrue_b32(), lo_sve, svext_u32(z, hi_sve, 4));
    }

    template <int n>
    static FORCE_INLINE XV shr128(const XV& a)
    {
        if constexpr (n == 0) return a;
        else if constexpr (n >= 16) return zero();
        svuint32_t z = svdup_u32(0);
        uint32x4_t lo = svget_neonq_u32(a.m_v);
        uint32x4_t hi = svget_neonq_u32(svext_u32(a.m_v, z, 4));
        lo = vreinterpretq_u32_u8(vextq_u8(vreinterpretq_u8_u32(lo), vdupq_n_u8(0), n));
        hi = vreinterpretq_u32_u8(vextq_u8(vreinterpretq_u8_u32(hi), vdupq_n_u8(0), n));
        svuint32_t lo_sve = svset_neonq_u32(z, lo);
        svuint32_t hi_sve = svset_neonq_u32(z, hi);
        return svorr_u32_x(svptrue_b32(), lo_sve, svext_u32(z, hi_sve, 4));
    }

    void broadcastLo128()
    {
        svuint32_t z = svdup_u32(0);
        uint32x4_t lo = svget_neonq_u32(m_v);
        svuint32_t lo_sve = svset_neonq_u32(z, lo);
        m_v = svorr_u32_x(svptrue_b32(), lo_sve, svext_u32(z, lo_sve, 4));
    }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b). Uses svbsl if SVE2 available.
    static FORCE_INLINE XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
    #if defined(__ARM_FEATURE_SVE2)
        return svbsl_u32_x(svptrue_b32(), a.m_v, b.m_v, mask.m_v);
    #else
        svbool_t p = svptrue_b32();
        return svorr_u32_x(p, svand_u32_x(p, mask.m_v, a.m_v), svbic_u32_x(p, b.m_v, mask.m_v));
    #endif
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
        return a ^ (b & c);
    }


    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        // sign-extend bit 0 of each 64-bit lane to a full-lane mask
        svuint64_t shifted = svlsl_n_u64_x(svptrue_b64(), svreinterpret_u64_u32(cond.m_v), 63);
        svint64_t  mask    = svasr_n_s64_x(svptrue_b64(), svreinterpret_s64_u64(shifted), 63);
        return sveor_u32_x(svptrue_b32(), m_v,
                           svand_u32_x(svptrue_b32(), svreinterpret_u32_s64(mask), cst.m_v));
    }

    // Conditional masking: For each 32-bit lane i, result[i] = (this[i] & 1) ? cst32[i] : 0.
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        svuint32_t shifted = svlsl_n_u32_x(svptrue_b32(), m_v, 31);
        svint32_t  mask    = svasr_n_s32_x(svptrue_b32(), svreinterpret_s32_u32(shifted), 31);
        return svand_u32_x(svptrue_b32(), svreinterpret_u32_s32(mask), cst32.m_v);
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    // Conditional XOR: For each 32-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    // Returns the parity (reduction XOR) of all bits in the register.
    uint8_t parity() const
    {
        svuint32_t z = svdup_u32(0);
        uint32x4_t lo = svget_neonq_u32(m_v);
        uint32x4_t hi = svget_neonq_u32(svext_u32(m_v, z, 4));
        return SimdRegister<128, ISA::NEON>(veorq_u32(lo, hi)).parity();
    }
};

#else  // x86 AVX2
template <>
struct SimdRegister<256, ISA::AVX2, void> : VirtualRegBase<256, ISA::AVX2>
{
    __m256i m_v;

    typedef SimdRegister<256, ISA::AVX2> XV;

    SimdRegister() : m_v(_mm256_undefined_si256()) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm256_load_si256((const __m256i*)p)) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm256_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(_mm256_set1_epi64x((long long)v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3)
    {
        __m128i tmp = _mm_setr_epi32(v0, v1, v2, v3);
        m_v = _mm256_set_m128i(tmp, tmp);
    }
    SimdRegister(__m256i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { if (A) _mm256_store_si256((__m256i*)dst, m_v); else _mm256_storeu_si256((__m256i*)dst, m_v); }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { if constexpr (A) return _mm256_load_si256((const __m256i*)p); else return _mm256_loadu_si256((const __m256i*)p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm256_or_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator*(const XV& a, const XV& b) { return _mm256_mullo_epi32(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm256_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm256_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm256_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm256_srli_epi64(a.m_v, n); }

    // Carry-less multiplication of two 64-bit polynomials (PCLMULQDQ).
    // imm8: 0x00: lo-lo, 0x01: hi-lo, 0x10: lo-hi, 0x11: hi-hi
    // Note: This operation is performed on each 128-bit lane.
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
#ifdef __VPCLMULQDQ__
        return _mm256_clmulepi64_si256(a.m_v, b.m_v, imm8);
#else
        // Emulate using 128-bit PCLMULQDQ
        __m128i a_lo = _mm256_castsi256_si128(a.m_v);
        __m128i a_hi = _mm256_extracti128_si256(a.m_v, 1);
        __m128i b_lo = _mm256_castsi256_si128(b.m_v);
        __m128i b_hi = _mm256_extracti128_si256(b.m_v, 1);
        __m128i res_lo = _mm_clmulepi64_si128(a_lo, b_lo, imm8);
        __m128i res_hi = _mm_clmulepi64_si128(a_hi, b_hi, imm8);
        return _mm256_set_m128i(res_hi, res_lo);
#endif
    }

    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
#ifdef __AVX512VL__
        __mmask8 isOdd = _mm256_test_epi64_mask(cond.m_v, _mm256_set1_epi64x(1LL));
        return _mm256_mask_xor_epi64(m_v, isOdd, m_v, cst.m_v);
#else
        __m256i lowestBit = _mm256_slli_epi64(cond.m_v, 63);
        __m256i hi   = _mm256_shuffle_epi32(lowestBit, 0xF5);
        __m256i mask = _mm256_srai_epi32(hi, 31);
        return _mm256_xor_si256(m_v, _mm256_and_si256(mask, cst.m_v));
#endif
    }

    // Concatenates registers a and b, then extracts a register-sized window starting from the n32-th word.
    // Effectively shifts the combined [a, b] window left by n32 words.
    // Example (n32FromSecond=1, 256-bit):
    //   a = {a0, a1, a2, a3, a4, a5, a6, a7}
    //   b = {b0, b1, b2, b3, b4, b5, b6, b7}
    //   result = {a1, a2, a3, a4, a5, a6, a7, b0}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 8, "n32FromSecond must be <= 8 for AVX2");

        if constexpr (n32FromSecond == 0)
            return a;
        else if constexpr (n32FromSecond == 8)
            return b;
#ifdef __AVX512VL__
        else {
            if constexpr (n32FromSecond == 1) return _mm256_alignr_epi32(b.m_v, a.m_v, 1);
            else if constexpr (n32FromSecond == 2) return _mm256_alignr_epi32(b.m_v, a.m_v, 2);
            else if constexpr (n32FromSecond == 3) return _mm256_alignr_epi32(b.m_v, a.m_v, 3);
            else if constexpr (n32FromSecond == 4) return _mm256_alignr_epi32(b.m_v, a.m_v, 4);
            else if constexpr (n32FromSecond == 5) return _mm256_alignr_epi32(b.m_v, a.m_v, 5);
            else if constexpr (n32FromSecond == 6) return _mm256_alignr_epi32(b.m_v, a.m_v, 6);
            else if constexpr (n32FromSecond == 7) return _mm256_alignr_epi32(b.m_v, a.m_v, 7);
            else return a; // should not happen due to if/else above
        }
#else
        else {
            // Combine the high 128 bits of a with the low 128 bits of b.
            __m256i aHibLo = _mm256_permute2x128_si256(a.m_v, b.m_v, 0x21);

            if constexpr (n32FromSecond < 4) {
                // Align: take bytes from a and then from combined.
                return _mm256_alignr_epi8(aHibLo, a.m_v, 4*n32FromSecond);
            }
            else if constexpr (n32FromSecond == 4) {
                return aHibLo;
            }
            else {
                // Need (n32FromSecond - 4) 32-bit words into the concatenation (b || aHibLo)
                return _mm256_alignr_epi8(b.m_v, aHibLo, 4*n32FromSecond - 16);
            }
        }
#endif
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm256_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm256_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm256_broadcastsi128_si256(_mm256_castsi256_si128(m_v)); }

    static FORCE_INLINE XV zero() { return _mm256_setzero_si256(); }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b). Uses ternary logic if AVX-512VL available.
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
#ifdef __AVX512VL__
        return _mm256_ternarylogic_epi32(mask.m_v, a.m_v, b.m_v, 0xCA);
#else
        return b ^ (mask & (a ^ b));
#endif
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
#ifdef __AVX512VL__
        return _mm256_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, 0x78);
#else
        return a ^ (b & c);
#endif
    }

#ifdef __AVX512VL__
    // General ternary logic operation (AVX-512VL).
    template <int imm>
    static FORCE_INLINE XV ternary(const XV& a, const XV& b, const XV& c)
    {
        return _mm256_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, imm);
    }
#endif

    // Conditional masking: For each 32-bit lane i, result[i] = (this[i] & 1) ? cst32[i] : 0.
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
#ifdef __AVX512VL__
        __mmask8 isOdd = _mm256_test_epi32_mask(m_v, _mm256_set1_epi32(1));
        return _mm256_maskz_mov_epi32(isOdd, cst32.m_v);
#else
        const __m256i z = _mm256_setzero_si256();
        const __m256i lowestBit = _mm256_slli_epi32(m_v, 31); // move least significant bit to most significant bit
        const __m256 mask = _mm256_castsi256_ps(lowestBit);
        return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(z), _mm256_castsi256_ps(cst32.m_v), mask));
#endif
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    // Conditional XOR: For each 32-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    // Returns the parity (reduction XOR) of all bits in the register.
    uint8_t parity() const
    {
        __m128i hi(_mm256_extracti128_si256(m_v, 1));
        __m128i lo(_mm256_castsi256_si128(m_v));
        return SimdRegister<128, ISA::SSE2>(_mm_xor_si128(lo, hi)).parity();
    }
};
#endif  // ARM_FEATURE_SVE vs x86 AVX2

#endif  // SIMD_N_BITS >= 256

#if SIMD_N_BITS>=512
template <>
struct SimdRegister<512, ISA::AVX512, void> : VirtualRegBase<512, ISA::AVX512>
{
    __m512i m_v;

    typedef SimdRegister<512, ISA::AVX512> XV;

    SimdRegister() : m_v(_mm512_setzero_si512()) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm512_load_si512((const __m512i*)p)) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm512_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(_mm512_set1_epi64((long long)v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm512_setr4_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(__m512i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) const { if (A) _mm512_store_si512((__m512i*)dst, m_v); else _mm512_storeu_si512((__m512i*)dst, m_v); }

    template <bool A = true>
    static FORCE_INLINE XV load(const void* p) { if constexpr (A) return _mm512_load_si512(p); else return _mm512_loadu_si512(p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm512_or_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm512_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm512_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm512_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm512_srli_epi64(a.m_v, n); }

    // Carry-less multiplication of two 64-bit polynomials (PCLMULQDQ).
    // imm8: 0x00: lo-lo, 0x01: hi-lo, 0x10: lo-hi, 0x11: hi-hi
    // Note: This operation is performed on each 128-bit lane.
    template <int imm8>
    FORCE_INLINE static XV clmul(const XV& a, const XV& b)
    {
#ifdef __VPCLMULQDQ__
        return _mm512_clmulepi64_si512(a.m_v, b.m_v, imm8);
#else
        // Emulate using 128-bit PCLMULQDQ
        __m128i a0 = _mm512_extracti32x4_epi32(a.m_v, 0);
        __m128i a1 = _mm512_extracti32x4_epi32(a.m_v, 1);
        __m128i a2 = _mm512_extracti32x4_epi32(a.m_v, 2);
        __m128i a3 = _mm512_extracti32x4_epi32(a.m_v, 3);
        __m128i b0 = _mm512_extracti32x4_epi32(b.m_v, 0);
        __m128i b1 = _mm512_extracti32x4_epi32(b.m_v, 1);
        __m128i b2 = _mm512_extracti32x4_epi32(b.m_v, 2);
        __m128i b3 = _mm512_extracti32x4_epi32(b.m_v, 3);
        __m128i r0 = _mm_clmulepi64_si128(a0, b0, imm8);
        __m128i r1 = _mm_clmulepi64_si128(a1, b1, imm8);
        __m128i r2 = _mm_clmulepi64_si128(a2, b2, imm8);
        __m128i r3 = _mm_clmulepi64_si128(a3, b3, imm8);
        __m512i res = _mm512_castsi128_si512(r0);
        res = _mm512_inserti32x4_epi32(res, r1, 1);
        res = _mm512_inserti32x4_epi32(res, r2, 2);
        res = _mm512_inserti32x4_epi32(res, r3, 3);
        return res;
#endif
    }

    // Conditional XOR: For each 64-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        __mmask8 isOdd = _mm512_test_epi64_mask(cond.m_v, _mm512_set1_epi64(1LL));
        return _mm512_mask_xor_epi64(m_v, isOdd, m_v, cst.m_v);
    }

    // Concatenates registers a and b, then extracts a register-sized window starting from the n32-th word.
    // Effectively shifts the combined [a, b] window left by n32 words.
    // Example (n32FromSecond=1, 512-bit):
    //   a = {a0, ..., a15}
    //   b = {b0, ..., b15}
    //   result = {a1, ..., a15, b0}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(XV a, XV b)
    {
        static_assert(n32FromSecond <= 16, "nBytesFromSecond must be <= 16 for AVX512");

        if constexpr (n32FromSecond == 0)
            return a;
        if constexpr (n32FromSecond == 16)
            return b;
        else {
            if constexpr (n32FromSecond == 1) return _mm512_alignr_epi32(b.m_v, a.m_v, 1);
            else if constexpr (n32FromSecond == 2) return _mm512_alignr_epi32(b.m_v, a.m_v, 2);
            else if constexpr (n32FromSecond == 3) return _mm512_alignr_epi32(b.m_v, a.m_v, 3);
            else if constexpr (n32FromSecond == 4) return _mm512_alignr_epi32(b.m_v, a.m_v, 4);
            else if constexpr (n32FromSecond == 5) return _mm512_alignr_epi32(b.m_v, a.m_v, 5);
            else if constexpr (n32FromSecond == 6) return _mm512_alignr_epi32(b.m_v, a.m_v, 6);
            else if constexpr (n32FromSecond == 7) return _mm512_alignr_epi32(b.m_v, a.m_v, 7);
            else if constexpr (n32FromSecond == 8) return _mm512_alignr_epi32(b.m_v, a.m_v, 8);
            else if constexpr (n32FromSecond == 9) return _mm512_alignr_epi32(b.m_v, a.m_v, 9);
            else if constexpr (n32FromSecond == 10) return _mm512_alignr_epi32(b.m_v, a.m_v, 10);
            else if constexpr (n32FromSecond == 11) return _mm512_alignr_epi32(b.m_v, a.m_v, 11);
            else if constexpr (n32FromSecond == 12) return _mm512_alignr_epi32(b.m_v, a.m_v, 12);
            else if constexpr (n32FromSecond == 13) return _mm512_alignr_epi32(b.m_v, a.m_v, 13);
            else if constexpr (n32FromSecond == 14) return _mm512_alignr_epi32(b.m_v, a.m_v, 14);
            else if constexpr (n32FromSecond == 15) return _mm512_alignr_epi32(b.m_v, a.m_v, 15);
            else return a;
        }
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm512_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm512_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm512_broadcast_i32x4(_mm512_castsi512_si128(m_v)); }

    // Conditional masking: For each 32-bit lane i, result[i] = (this[i] & 1) ? cst32[i] : 0.
    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        const __mmask16 isOdd = _mm512_test_epi32_mask(m_v, _mm512_set1_epi32(1));
        return _mm512_maskz_mov_epi32(isOdd, cst32.m_v);
    }

    FORCE_INLINE XV ifOddCstThenZero(const XV cst) const { return ifOddCst32ElseZero(cst); }

    // General ternary logic operation (AVX-512).
    template <int imm>
    static FORCE_INLINE XV ternary(const XV& a, const XV& b, const XV& c)
    {
        return _mm512_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, imm);
    }

    static FORCE_INLINE XV zero() { return _mm512_setzero_si512(); }

    // Bitwise Selection: For each bit, result = (mask & a) | (~mask & b). Uses ternary logic.
    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return _mm512_ternarylogic_epi32(mask.m_v, a.m_v, b.m_v, 0xCA);
    }

    // result = a ^ (b & c)
    FORCE_INLINE static XV bitwiseXorAnd(const XV a, const XV b, const XV c)
    {
        return _mm512_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, 0x78);
    }

    // Conditional XOR: For each 32-bit lane i, result[i] = m_v[i] ^ (cond[i] & 1 ? cst[i] : 0).
    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        __mmask16 isOdd = _mm512_test_epi32_mask(cond.m_v, _mm512_set1_epi32(1));
        return _mm512_mask_xor_epi32(m_v, isOdd, m_v, cst.m_v);
    }

    // Returns the parity (reduction XOR) of all bits in the register.
    uint8_t parity() const
    {
        __m256i hi(_mm512_extracti64x4_epi64(m_v, 1));
        __m256i lo(_mm512_castsi512_si256(m_v));
        return SimdRegister<256, ISA::AVX2>(_mm256_xor_si256(lo, hi)).parity();
    }
};
#endif

} // namespace details
} // namespace xvmt
