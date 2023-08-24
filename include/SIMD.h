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
    FORCE_INLINE void store(uint32_t* p)
    {
        for (size_t i = 0; i < s_M; ++i, p += sizeof(XVHw) / sizeof(uint32_t))
            m_v[i].template store<A>(p);
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

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = m_v[i].xorIfOddCst64(cond.m_v[i], cst.m_v[i]);
        return r;
    }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = XVHw::bitwiseSelect(mask.m_v[i], a.m_v[i], b.m_v[i]);
        return r;
    }

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
                XVHw current = a.m_v[4 * s];
                for (size_t j = 0; j < 3; ++j) {
                    XVHw next = a.m_v[4 * s +  j + 1];
                    r.m_v[4 * s + j] = (current >> 8) | (next << 24);
                    current = next;
                }
                r.m_v[4 * s + 3] = current >> 8;
            }
        }
        else {
            for (size_t i = 0; i < s_M; ++i)
                r.m_v[i] = XVHw::template shr128<nBytes>(a.m_v[i]);
        }
        return r;
    }

    FORCE_INLINE void broadcastLo128()
    {
        static_assert(N128 > 0);
        if constexpr (HwBitLen == 32) {
            for (size_t i = 0; i < N128; ++i)
                for (size_t j = 0; j < 4; ++j)
                    m_v[4 * i + j] = m_v[j];
        }
        else {
            for (size_t i = 1; i < s_M; ++i)
                m_v[i] = m_v[0];
        }
    }

    bool eq(const XV& rhs) const
    {
        for (size_t i = 0; i < s_M; ++i)
            if (m_v[i] != rhs.m_v[i])
                return false;
        return true;
    }

    FORCE_INLINE static XV zero()
    {
        XVHw z(0);
        XV r;
        for (auto& v : r.m_v)
            v = z;
        return r;
    }

    template <typename XVI>
    FORCE_INLINE XV ifOddCst32ElseZero(const XVI value) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i) {
            if constexpr (XVI::s_virtualBitLen > HwBitLen)
                r.m_v[i] = m_v[i].ifOddCst32ElseZero(value.m_v[i]);
            else
                r.m_v[i] = m_v[i].ifOddCst32ElseZero(value);
        }
        return r;
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        XV r;
        for (size_t i = 0; i < s_M; ++i)
            r.m_v[i] = m_v[i].xorIfOddCst32(cond.m_v[i], cst.m_v[i]);
        return r;
    }

    FORCE_INLINE uint8_t parity() const
    {
        XVHw temp(m_v[0]);
        for (size_t i = 1; i < s_M; ++i)
            temp = temp ^ m_v[i];
        return temp.parity();
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
    FORCE_INLINE void store(uint32_t* dst) { *dst = m_v; }

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return a.m_v > b.m_v ? uint32_t(-1) : uint32_t(0); }
    friend FORCE_INLINE XV operator>>(const XV a, int n) { return uint32_t(a.m_v >> n); }
    friend FORCE_INLINE XV operator<<(const XV a, int n) { return uint32_t(a.m_v << n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
#if 0 && defined(__GNUC__) && defined(__x86_64__)
        // force the use of cmov with gcc
        uint32_t z;
        __asm__(
            "mov %[a], %[z]\n"
            "and $0x1, %[z]\n"
            "cmovne %[b], %[z]\n"
            : [z] "=r"(z)
            : [a] "r"(m_v), [b] "r"(cst32.m_v)
            : "cc"
        );
        return z;
#elif defined(_MSC_VER)
        const uint32_t lowestBit = m_v & 0x1;
        return lowestBit ? cst32 : zero();
#else
        const uint32_t x[2] = { 0, cst32.m_v };
        const uint32_t lowestBit = m_v & 0x1;
        return x[lowestBit];
#endif
    }

    // combine: {x0, x1, x2, x3} -> {x4, x5, x6, x7} => {x1, x2, x3, x4}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(XV a, XV b)
    {
        static_assert(n32FromSecond <= 1, "n32FromSecond must be <=1 with 32-bit registers");
        if constexpr (n32FromSecond == 0)
            return a;
        else if constexpr (n32FromSecond == 1)
            return b;
    }

    static FORCE_INLINE XV zero() { return uint32_t(0); }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return (mask.m_v & a.m_v) | (~mask.m_v & b.m_v);
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
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
    FORCE_INLINE void store(uint32_t* dst) { *(uint64_t*)dst = m_v; }

    friend FORCE_INLINE XV operator&(const XV a, const XV b) { return a.m_v & b.m_v; }
    friend FORCE_INLINE XV operator^(const XV a, const XV b) { return a.m_v ^ b.m_v; }
    friend FORCE_INLINE XV operator|(const XV a, const XV b) { return a.m_v | b.m_v; }

    friend FORCE_INLINE XV shr64(const XV a, int n) { return uint64_t(a.m_v >> n); }
    friend FORCE_INLINE XV shl64(const XV a, int n) { return uint64_t(a.m_v << n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return m_v == rhs.m_v; }

    static FORCE_INLINE XV zero() { return uint64_t(0); }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return (mask.m_v & a.m_v) | (~mask.m_v & b.m_v);
    }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        return *this ^ XV(cond.m_v & 1 ? cst.m_v : uint64_t(0));
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
    FORCE_INLINE void store(uint32_t* dst) { vst1q_u32(dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return vandq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return veorq_u32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return vorrq_u32(a.m_v, b.m_v); }
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

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return vbslq_u32(mask.m_v, a.m_v, b.m_v);
    }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        int32x4_t mask = vshrq_n_s32(vreinterpretq_s32_u32(vshlq_n_u32(m_v, 31)), 31);
        return vandq_u32(vreinterpretq_u32_s32(mask), cst32.m_v);
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        // sign-extend bit 0 of each 64-bit lane to a full-lane mask
        int64x2_t mask = vshrq_n_s64(vreinterpretq_s64_u64(vshlq_n_u64(vreinterpretq_u64_u32(cond.m_v), 63)), 63);
        return veorq_u32(m_v, vandq_u32(vreinterpretq_u32_s64(mask), cst.m_v));
    }

    uint8_t parity() const
    {
        uint32_t d = vgetq_lane_u32(m_v, 0) ^ vgetq_lane_u32(m_v, 1) ^
                     vgetq_lane_u32(m_v, 2) ^ vgetq_lane_u32(m_v, 3);
        return popcnt(d) & 1;
    }
};
#else
template <ISA Isa>
struct SimdRegister<128, Isa, std::enable_if_t<Isa == ISA::SSE2 || Isa == ISA::SSE42, void>> : VirtualRegBase<128, Isa>
{
    __m128i m_v;

    typedef SimdRegister<128, Isa> XV;

    SimdRegister() : m_v(_mm_undefined_si128()) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(_mm_set1_epi64x((long long)v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm_setr_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm_load_si128((const __m128i*)p)) {}
    FORCE_INLINE SimdRegister(__m128i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm_store_si128((__m128i*)dst, m_v); else _mm_storeu_si128((__m128i*)dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm_and_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm_xor_si128(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm_or_si128(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm_srli_epi64(a.m_v, n); }

    FORCE_INLINE bool eq(const XV& rhs) const { return _mm_test_all_ones(_mm_cmpeq_epi32(m_v, rhs.m_v)); }

    // combine: {x0, x1, x2, x3} -> {x4, x5, x6, x7} => {x1, x2, x3, x4}
    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        // Use SSSE3/_mm_alignr_epi8 to concatenate (b | a) and extract starting at byte nBytesFromSecond
        // This drops the lowest nBytesFromSecond bytes from a and append the nBytesFromSecond bytes from b
        // For example, if nBytesFromSecond = 4
        // combine<4>: {x0, x1, x2, x3} -> {x4, x5, x6, x7} => {x1, x2, x3, x4}
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

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return _mm_or_si128(_mm_and_si128(mask.m_v, a.m_v), _mm_andnot_si128(mask.m_v, b.m_v));
    }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
#if 0
        const __m128 z = _mm_setzero_ps();
        const __m128 lowestBit = _mm_castsi128_ps(_mm_slli_epi32(m_v, 31));
        return _mm_castps_si128(_mm_blendv_ps(z, _mm_castsi128_ps(cst32.m_v), lowestBit));
#else
        const __m128i z = zero().m_v;
        const __m128i lowestBit = _mm_slli_epi32(m_v, 31);
        const __m128i isOdd = _mm_cmpgt_epi32(z, lowestBit);
        return _mm_and_si128(isOdd, cst32.m_v);
#endif
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        // sign-extend bit 0 of each 64-bit lane to a full-lane mask
        __m128i lowestBit = _mm_slli_epi64(cond.m_v, 63);
        __m128i hi   = _mm_shuffle_epi32(lowestBit, 0xF5);  // broadcast high-32 of each 64-bit lane
        __m128i mask = _mm_srai_epi32(hi, 31);
        return _mm_xor_si128(m_v, _mm_and_si128(mask, cst.m_v));
    }

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
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm256_store_si256((__m256i*)dst, m_v); else _mm256_storeu_si256((__m256i*)dst, m_v); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm256_and_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm256_xor_si256(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm256_or_si256(a.m_v, b.m_v); }
    //friend FORCE_INLINE XV operator>(const XV& a, const XV& b) { return _mm256_cmpgt_epi32(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm256_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm256_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm256_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm256_srli_epi64(a.m_v, n); }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        __m256i lowestBit = _mm256_slli_epi64(cond.m_v, 63);
        __m256i hi   = _mm256_shuffle_epi32(lowestBit, 0xF5);
        __m256i mask = _mm256_srai_epi32(hi, 31);
        return _mm256_xor_si256(m_v, _mm256_and_si256(mask, cst.m_v));
    }

    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(const XV& a, const XV& b)
    {
        static_assert(n32FromSecond <= 8, "n32FromSecond must be <= 8 for AVX2");

        if constexpr (n32FromSecond == 0)
            return a;
        else if constexpr (n32FromSecond == 8)
            return b;
        else {
            // Combine the high 128 bits of v0 with the low 128 bits of v1.
            __m256i aHibLo = _mm256_permute2x128_si256(a.m_v, b.m_v, 0x21);

            if constexpr (n32FromSecond < 4) {
                // Align: take bytes from v0 and then from combined.
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
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm256_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm256_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm256_broadcastsi128_si256(_mm256_castsi256_si128(m_v)); }

    static FORCE_INLINE XV zero() { return _mm256_setzero_si256(); }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(b.m_v), _mm256_castsi256_ps(a.m_v), _mm256_castsi256_ps(mask.m_v)));
    }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        const __m256i z = zero().m_v;
        const __m256i lowestBit = _mm256_slli_epi32(m_v, 31); // move least significant bit to most significant bit
#if 0
        const __m256i isOdd = _mm256_cmpgt_epi32(z, lowestBit);
        return _mm256_and_si256(value.m_v, isOdd);
#else
        const __m256 mask = _mm256_castsi256_ps(lowestBit);
        return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(z), _mm256_castsi256_ps(cst32.m_v), mask));
#endif
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        return *this ^ cond.ifOddCst32ElseZero(cst);
    }

    uint8_t parity() const
    {
        __m128i hi(_mm256_extracti128_si256(m_v, 1));
        __m128i lo(_mm256_castsi256_si128(m_v));
        return SimdRegister<128, ISA::SSE2>(_mm_xor_si128(lo, hi)).parity();
    }
};
#endif

#if SIMD_N_BITS>=512
template <>
struct SimdRegister<512, ISA::AVX512, void> : VirtualRegBase<512, ISA::AVX512>
{
    __m512i m_v;

    typedef SimdRegister<512, ISA::AVX512> XV;

    SimdRegister() : m_v(_mm512_undefined_si512()) {}
    FORCE_INLINE SimdRegister(const void* p) : m_v(_mm512_load_si512((const __m512i*)p)) {}
    FORCE_INLINE SimdRegister(uint32_t v) : m_v(_mm512_set1_epi32(v)) {}
    FORCE_INLINE SimdRegister(uint64_t v) : m_v(_mm512_set1_epi64((long long)v)) {}
    FORCE_INLINE SimdRegister(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) : m_v(_mm512_setr4_epi32(v0, v1, v2, v3)) {}
    FORCE_INLINE SimdRegister(__m512i v) : m_v(v) {}

    template <bool A>
    FORCE_INLINE void store(uint32_t* dst) { if (A) _mm512_store_si512((__m512i*)dst, m_v); else _mm512_storeu_si512((__m512i*)dst, m_v); }

    //    static FORCE_INLINE XV load(const void* p) { return _mm512_loadu_si512((const __m512i*) p); }

    friend FORCE_INLINE XV operator&(const XV& a, const XV& b) { return _mm512_and_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator^(const XV& a, const XV& b) { return _mm512_xor_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator|(const XV& a, const XV& b) { return _mm512_or_si512(a.m_v, b.m_v); }
    friend FORCE_INLINE XV operator<<(const XV& a, const int n) { return _mm512_slli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV operator>>(const XV& a, const int n) { return _mm512_srli_epi32(a.m_v, n); }
    friend FORCE_INLINE XV shl64(const XV& a, int n) { return _mm512_slli_epi64(a.m_v, n); }
    friend FORCE_INLINE XV shr64(const XV& a, int n) { return _mm512_srli_epi64(a.m_v, n); }

    FORCE_INLINE XV xorIfOddCst64(const XV& cond, const XV& cst) const
    {
        __mmask8 isOdd = _mm512_test_epi64_mask(cond.m_v, _mm512_set1_epi64(1LL));
        return _mm512_mask_xor_epi64(m_v, isOdd, m_v, cst.m_v);
    }

    template <unsigned n32FromSecond>
    static FORCE_INLINE XV alignr32(XV a, XV b)
    {
        static_assert(n32FromSecond <= 16, "nBytesFromSecond must be <= 16 for AVX512");

        if constexpr (n32FromSecond == 0)
            return a;
        if constexpr (n32FromSecond == 16)
            return b;
        else {
            // Combine the high 256 bits of 'a' with the low 256 bits of 'b'
            return _mm512_alignr_epi32(b.m_v, a.m_v, n32FromSecond);
        }
    }

    template <int n>
    static FORCE_INLINE XV shl128(const XV& a) { return _mm512_bslli_epi128(a.m_v, n); }
    template <int n>
    static FORCE_INLINE XV shr128(const XV& a) { return _mm512_bsrli_epi128(a.m_v, n); }

    void broadcastLo128() { m_v = _mm512_broadcast_i32x4(_mm512_castsi512_si128(m_v)); }

    FORCE_INLINE XV ifOddCst32ElseZero(const XV cst32) const
    {
        const __mmask16 isOdd = _mm512_test_epi32_mask(m_v, _mm512_set1_epi32(1));
        return _mm512_maskz_mov_epi32(isOdd, cst32.m_v);
    }

    static FORCE_INLINE XV ternary(const XV& a, const XV& b, const XV& c, int imm)
    {
        return _mm512_ternarylogic_epi32(a.m_v, b.m_v, c.m_v, imm);
    }

    static FORCE_INLINE XV zero() { return _mm512_setzero_si512(); }

    FORCE_INLINE static XV bitwiseSelect(const XV mask, const XV a, const XV b)
    {
        return _mm512_ternarylogic_epi32(mask.m_v, a.m_v, b.m_v, 0xCA);
    }

    FORCE_INLINE XV xorIfOddCst32(const XV& cond, const XV& cst) const
    {
        __mmask16 isOdd = _mm512_test_epi32_mask(cond.m_v, _mm512_set1_epi32(1));
        return _mm512_mask_xor_epi32(m_v, isOdd, m_v, cst.m_v);
    }

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

