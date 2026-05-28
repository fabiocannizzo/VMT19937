#ifndef XVMT_JUMP_TRAITS_H
#define XVMT_JUMP_TRAITS_H

#include "RandGen.h"
#include "polynomial_jump.h"
#include "jump_matrix.h"
#include "../SFMT-src-1.5.1/SFMT.h"

#include <string>
#include <filesystem>
#include <iostream>

extern "C" {
    void               init_genrand(unsigned long s);
    unsigned long      genrand_int32();
    void               init_genrand64(unsigned long long seed);
    unsigned long long genrand64_int64();
}

namespace xvmt {
namespace details {

struct MT32Traits {
    using Gen        = XMT19937<SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = MT19937Matrix<32>;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = XMT19937<ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint32_t;

    static constexpr int  s_reduceDegree  = 19937;
    static constexpr int  s_startStep     = 0;
    static constexpr size_t s_power2      = 5; // 2^5 = 32
    static constexpr const char* s_gentype       = "mt32";
    static constexpr const char* s_chainPrefix   = "mt32_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/mt32/characteristic.mt32.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/mt32/";

    static constexpr uint32_t s_seed = 1234UL;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/mt32/J%05d.mt32.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/mt32/F%05d.mt32.bits", step); return buf;
    }

    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seed, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seed, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seed, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seed, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint32(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint32(); }
};

struct MT64Traits {
    using Gen        = XMT19937_64<SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = MT19937Matrix<64>;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = XMT19937_64<ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint64_t;

    static constexpr int  s_reduceDegree  = 19937;
    static constexpr int  s_startStep     = 0;
    static constexpr size_t s_power2      = 6; // 2^6 = 64
    static constexpr const char* s_gentype       = "mt64";
    static constexpr const char* s_chainPrefix   = "mt64_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/mt64/characteristic.mt64.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/mt64/";

    static constexpr uint64_t s_seed = 0x123456789ABCULL;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/mt64/J%05d.mt64.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/mt64/F%05d.mt64.bits", step); return buf;
    }

    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seed, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seed, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seed, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seed, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint64(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint64(); }
};

struct SFMTTraits {
    using Gen        = VSFMT19937<128, false, SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = SFMT19937Matrix;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = VSFMT19937<128, false, ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint32_t;

    static constexpr int  s_reduceDegree  = 19968;
    static constexpr int  s_startStep     = 2;
    static constexpr size_t s_power2      = 2; // SFMT has 4 output words per recurrence step
    static constexpr const char* s_gentype       = "sfmt";
    static constexpr const char* s_chainPrefix   = "sfmt_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/sfmt/characteristic.sfmt.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/sfmt/";

    static constexpr uint32_t s_seeds[]  = { 0x123, 0x234, 0x345, 0x456 };
    static constexpr uint32_t s_seedLen  = 4;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/sfmt/J%05d.sfmt.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/sfmt/F%05d.sfmt.bits", step); return buf;
    }

    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seeds, s_seedLen, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seeds, s_seedLen, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seeds, s_seedLen, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seeds, s_seedLen, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint32(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint32(); }
};

} // namespace details
} // namespace xvmt

#endif // XVMT_JUMP_TRAITS_H
