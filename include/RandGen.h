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
    using output_word_t = typename base_t::output_word_t;

    // re-export useful constants
    static constexpr size_t s_regLenBits = base_t::s_regLenBits;           // logical SIMD register width in bits
    static constexpr size_t s_regLenBitsHw = base_t::s_regLenBitsHw;      // hardware SIMD register width in bits
    static constexpr size_t s_n32InOneWord = base_t::s_n32InOneWord;       // uint32 elements per MT/SFMT word
    static constexpr size_t s_n32InFullState = base_t::s_n32InFullState;   // total uint32 elements in the interleaved state array
    static constexpr size_t s_nStates = base_t::s_nStates;                 // number of parallel generator states

private:
    void completeStateInitialization(size_t nCommonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        // temporary workspace matrix
        BinaryMatrix<2, base_t::s_nMatrixBits> tmp(commonJump || sequentialJump);

        // apply common jump to state-0
        if (nCommonJumpRepeat) {
            MYASSERT(commonJump, "commonJump is required when nCommonJumpRepeat>0");

            // extract state 0
            base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

            for (size_t i = 0; i < nCommonJumpRepeat; ++i)
                commonJump->multiplyByColumn(tmp.rowBegin((i + 1) % 2), tmp.rowBegin(i % 2));

            // copy back to state 0
            base_t::vectorToState(0, (const uint32_t*)tmp.rowBegin(nCommonJumpRepeat % 2));
        }

        // if there are multiple states, distance them using the sequentialJump matrix
        if constexpr (s_nStates > 1) {
            if (sequentialJump) {
                // perform jump ahead of the s_regLenWords states
                // State_0 = State_0
                // State_1 = Jump x State_0
                // State_2 = Jump x State_1
                // ...

                // copy state 0 to the first row
                base_t::stateToVector(0, (uint32_t*)tmp.rowBegin(0));

                for (size_t s = 1; s < base_t::s_nStates; ++s) {
                    const uint8_t* psrc = (uint8_t*)tmp.rowBegin((s + 1) % 2);
                    uint8_t* pdst = (uint8_t*)tmp.rowBegin(s % 2);

                    sequentialJump->multiplyByColumn(pdst, psrc);

                    // copy to state vector s
                    base_t::vectorToState(s, (const uint32_t*)pdst);
                }
            }
            else {
#if (RANDGEN_TESTING!=1)
                THROW("Having multiple states and no sequential jump matrix does not make sense");
#endif
                // Fallback for RANDGEN_TESTING only: replicate state 0 to all states.
                // 32-bit output_word_t (MT32, SFMT) uses s_n32inReg-based interleaving (j-loop needed for SFMT);
                // 64-bit output_word_t (MT64) uses a direct word-level copy.
                if constexpr (sizeof(output_word_t) == 4) {
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

    // Initialize as follows:
    // 1) initialize state 0 with seed
    // 2) apply commonJump matrix nCommonJumpRepeat times to state 0
    // 3) if multiple states are present, apply sequentialJump matrix to initialize the other states
    // Note that the sequentialJump must be provided only for generators of the V-family, which have multiple states
    RandGen(output_word_t seed, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seed, commonJumpRepeat, commonJump, sequentialJump);
    }

    // Initialize as follows:
    // 1) initialize state 0 with seeds array
    // 2) apply commonJump matrix nCommonJumpRepeat times to state 0
    // 3) if multiple states are present, apply sequentialJump matrix to initialize the other states
    // Note that the sequentialJump must be provided only for generators of the V-family, which have multiple states
    RandGen(const uint32_t seeds[], uint32_t n_seeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
        : base_t()
    {
        reinit(seeds, n_seeds, commonJumpRepeat, commonJump, sequentialJump);
    }

    // Re-initialize as follows:
    // 1) initialize state 0 with seed
    // 2) apply commonJump matrix nCommonJumpRepeat times to state 0
    // 3) if multiple states are present, apply sequentialJump matrix to initialize the other states
    // Note that the sequentialJump must be provided only for generators of the V-family, which have multiple states
    void reinit(output_word_t s, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(s);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    // Re-initialize as follows:
    // 1) initialize state 0 with seeds array
    // 2) apply commonJump matrix nCommonJumpRepeat times to state 0
    // 3) if multiple states are present, apply sequentialJump matrix to initialize the other states
    // Note that the sequentialJump must be provided only for generators of the V-family, which have multiple states
    void reinit(const uint32_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(seeds, nSeeds);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    // Re-initialize as follows:
    // 1) initialize state 0 with seeds array
    // 2) apply commonJump matrix nCommonJumpRepeat times to state 0
    // 3) if multiple states are present, apply sequentialJump matrix to initialize the other states
    // Note that the sequentialJump must be provided only for generators of the V-family, which have multiple states
    void reinit(const uint64_t* seeds, uint32_t nSeeds, size_t commonJumpRepeat, const matrix_t* commonJump, const matrix_t* sequentialJump)
    {
        base_t::reinitMainState(seeds, nSeeds);
        completeStateInitialization(commonJumpRepeat, commonJump, sequentialJump);
        base_t::reinitPointers();
    }

    // generates a random number on [0,0xffffffff] interval
    FORCE_INLINE uint32_t genrand_uint32()
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Scalar or QM_Any");
        return base_t::genrand_uint32();
    }

    // generates a random number on [0,0xffffffffffffffff] interval
    FORCE_INLINE uint64_t genrand_uint64()
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Scalar or QM_Any");
        return base_t::genrand_uint64();
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QryBlk16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_uint32_blk16(dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    FORCE_INLINE void genrand_word_blk(output_word_t* dst)
    {
        static_assert(QryBlk16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_word_blk(dst);
    }

    // generates n uniform discrete random numbers in [0,0xffffffff] interval
    FORCE_INLINE void genrand_uint32_anySize(uint32_t* dst, size_t n)
    {
        static_assert(!QryBlk16, "This function can only be invoked when query mode is QM_Any");
        base_t::genrand_uint32_anySize(dst, n);
    }

    // generates n uniform discrete output_word_t random numbers
    FORCE_INLINE void genrand_word_anySize(output_word_t* dst, size_t n)
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
    using base_t::RandGen; // reuse constructors
};

template < ISA Isa = details::BestIsa<512>::isa
         , bool QryBlk16 = false
         >
struct XMT19937 : details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<IsaTraits<Isa>::HwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::isa
         >
struct VMT19937_64 : details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

template < ISA Isa = details::BestIsa<512>::isa
         , bool QryBlk16 = false
         >
struct XMT19937_64 : details::RandGen<details::MT19937Base<(IsaTraits<Isa>::HwBitLen >= 64 ? IsaTraits<Isa>::HwBitLen : 64), Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<(IsaTraits<Isa>::HwBitLen >= 64 ? IsaTraits<Isa>::HwBitLen : 64), Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::isa
         >
struct VSFMT19937 : details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>
{
    using base_t = details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

} // namespace xvmt
