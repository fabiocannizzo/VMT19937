#pragma once

#include "MT19937.h"
#include "SFMT19937.h"

namespace xvmt {
namespace details {

template <typename GenBase, bool QryBlk16>
class RandGen : protected GenBase
{
    using base_t = GenBase;
public:
    using matrix_t = typename base_t::matrix_t;
    using word_t = typename base_t::word_t;

    // re-export useful constants
    static constexpr size_t s_regLenBits = base_t::s_regLenBits;
    static constexpr size_t s_regLenBitsHw = base_t::s_regLenBitsHw;
    static constexpr size_t s_n32InOneWord = base_t::s_n32InOneWord;
    static constexpr size_t s_n32InFullState = base_t::s_n32InFullState;
    static constexpr size_t s_nStates = base_t::s_nStates;

private:
    void completeStateInitialization(size_t nCommonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        BinaryMatrix<2, base_t::s_nMatrixBits> tmp(commonJump || sequentialJump);

        if (nCommonJumpRepeat) {
            MYASSERT(commonJump, "commonJump is required when nCommonJumpRepeat>0");

            base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

            for (size_t i = 0; i < nCommonJumpRepeat; ++i)
                commonJump->multiplyByColumn(tmp.rowBegin((i + 1) % 2), tmp.rowBegin(i % 2));

            base_t::vectorToState(0, (const uint32_t*)tmp.rowBegin(nCommonJumpRepeat % 2));
        }

        if constexpr (s_nStates > 1) {
            if (sequentialJump) {
                base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

                for (size_t s = 1; s < base_t::s_nStates; ++s) {
                    const uint8_t* psrc = (uint8_t*)tmp.rowBegin((s + 1) % 2);
                    uint8_t* pdst = (uint8_t*)tmp.rowBegin(s % 2);

                    sequentialJump->multiplyByColumn(pdst, psrc);

                    base_t::vectorToState(s, (const uint32_t*)pdst);
                }
            }
            else {
#if (RANDGEN_TESTING!=1)
                THROW("Having multiple states and no sequential jump matrix does not make sense");
#endif
                // Fallback for RANDGEN_TESTING only: replicate state 0 to all states.
                // 32-bit word_t (MT32, SFMT) uses s_n32inReg-based interleaving (j-loop needed for SFMT);
                // 64-bit word_t (MT64) uses a direct word-level copy.
                if constexpr (sizeof(word_t) == 4) {
                    for (size_t w = 0; w < (size_t)base_t::s_N; ++w)
                        for (size_t j = 0; j < base_t::s_n32InOneWord; ++j)
                            for (size_t s = 1; s < base_t::s_nStates; ++s)
                                base_t::m_state[w * base_t::s_n32inReg + s * base_t::s_n32InOneWord + j] = base_t::m_state[w * base_t::s_n32inReg + j];
                } else {
                    for (size_t w = 0; w < (size_t)base_t::s_N; ++w)
                        for (size_t s = 1; s < base_t::s_nStates; ++s)
                            base_t::m_state[w * base_t::s_nStates + s] = base_t::m_state[w * base_t::s_nStates];
                }
            }
        }
        else {
            MYASSERT(!sequentialJump, "sequentialJump matrix should not be provided when there is only one state");
        }
    }

public:
    RandGen() {}

    RandGen(word_t seed, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    RandGen(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seeds, n_seeds, commonJumpRepeat, commonJump, sequentialJump);
    }

    void reinit(word_t s, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(s);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    void reinit(const uint32_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(seeds, nSeeds);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    void reinit(const uint64_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(seeds, nSeeds);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Scalar or QM_Any");
        return base_t::genrand_uint32();
    }

    FORCE_INLINE uint64_t genrand_uint64()
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Scalar or QM_Any");
        return base_t::genrand_uint64();
    }

    FORCE_INLINE void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QryBlk16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_uint32_blk16(dst);
    }

    FORCE_INLINE void genrand_word_blk(word_t* dst)
    {
        static_assert(QryBlk16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_word_blk(dst);
    }

    FORCE_INLINE void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Any");
        base_t::genrand_uint32_anySize(dst, n);
    }

    FORCE_INLINE void genrand_word_anySize(word_t* dst, size_t n)
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Any");
        base_t::genrand_word_anySize(dst, n);
    }
};

} // namespace details

template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::isa
         >
struct VMT19937 : details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<32>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<32>>, QryBlk16>;
    using base_t::RandGen;
};

template < ISA Isa = details::BestIsa<512>::isa
         , bool QryBlk16 = false
         >
struct XMT19937 : details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>;
    using base_t::RandGen;
};

template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::isa
         >
struct VMT19937_64 : details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen;
};

template < ISA Isa = details::BestIsa<512>::isa
         , bool QryBlk16 = false
         >
struct XMT19937_64 : details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen;
};

template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::isa
         >
struct VSFMT19937 : details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>
{
    using base_t = details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>;
    using base_t::RandGen;
};

} // namespace xvmt
