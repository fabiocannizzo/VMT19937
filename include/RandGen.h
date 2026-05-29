#pragma once

#include "MT19937.h"
#include "SFMT19937.h"
#include "polynomial_jump.h"
#include <type_traits>

namespace xvmt {
namespace details {

template <typename GenBase, bool QryBlk16>
class RandGen : protected GenBase
{
    using base_t = GenBase;
public:
    using matrix_t = typename base_t::matrix_t;
    using output_word_t = typename base_t::output_word_t;
    using poly_t = Polynomial<32768, base_t::s_isa>;

    // re-export useful constants
    static constexpr size_t s_regLenBits = base_t::s_regLenBits;           // logical SIMD register width in bits
    static constexpr size_t s_regLenBitsHw = base_t::s_regLenBitsHw;      // hardware SIMD register width in bits
    static constexpr size_t s_n32InOneWord = base_t::s_n32InOneWord;       // uint32 elements per MT/SFMT word
    static constexpr size_t s_n32InFullState = base_t::s_n32InFullState;   // total uint32 elements in the interleaved state array
    static constexpr size_t s_nStates = base_t::s_nStates;                 // number of parallel generator states
    static constexpr BitsGenType s_bitsGenType = base_t::s_bitsGenType;

private:
    /**
     * @brief Policy for matrix-based state advancement.
     */
    struct MatrixJumpPolicy {
        size_t nRepeat;
        const matrix_t* common;
        const matrix_t* sequential;

        bool hasCommon() const { return nRepeat > 0 && common != nullptr; }
        bool hasSequential() const { return sequential != nullptr; }

        void applyCommon(RandGen& gen) const {
            BinaryMatrix<2, base_t::s_nMatrixBits> tmp(true);
            gen.stateToVector(0, (uint32_t*)tmp.rowBegin(0));
            for (size_t i = 0; i < nRepeat; ++i)
                common->multiplyByColumn(tmp.rowBegin((i + 1) % 2), tmp.rowBegin(i % 2));
            gen.vectorToState(0, (const uint32_t*)tmp.rowBegin(nRepeat % 2));
        }

        void applySequential(RandGen& gen) const {
            BinaryMatrix<2, base_t::s_nMatrixBits> tmp(true);
            gen.stateToVector(0, (uint32_t*)tmp.rowBegin(0));
            for (size_t s = 1; s < s_nStates; ++s) {
                const uint8_t* psrc = (uint8_t*)tmp.rowBegin((s + 1) % 2);
                uint8_t* pdst = (uint8_t*)tmp.rowBegin(s % 2);
                sequential->multiplyByColumn(pdst, psrc);
                gen.vectorToState(s, (const uint32_t*)pdst);
            }
        }
    };

    /**
     * @brief Policy for polynomial-based state advancement.
     */
    struct PolyJumpPolicy {
        const poly_t* common;
        const poly_t* sequential;

        bool hasCommon() const { return common != nullptr; }
        bool hasSequential() const { return sequential != nullptr; }

        void applyCommon(RandGen& gen) const {
            PolynomialJumpApplier<base_t>::apply(gen, *common, 0);
        }

        void applySequential(RandGen& gen) const {
            for (size_t s = 1; s < s_nStates; ++s) {
                gen.copyState(s - 1, s);
                PolynomialJumpApplier<base_t>::apply(gen, *sequential, s);
            }
        }
    };

    /**
     * @brief Copy the internal state of one parallel instance to another.
     */
    void copyState(size_t srcIdx, size_t dstIdx) {
        if constexpr (sizeof(output_word_t) == 4) {
            for (size_t w = 0; w < (size_t)base_t::s_N; ++w)
                for (size_t j = 0; j < base_t::s_n32InOneWord; ++j)
                    base_t::m_state[w * base_t::s_n32inReg + dstIdx * base_t::s_n32InOneWord + j] = 
                        base_t::m_state[w * base_t::s_n32inReg + srcIdx * base_t::s_n32InOneWord + j];
        } else {
            for (size_t w = 0; w < (size_t)base_t::s_N; ++w)
                base_t::m_state[w * base_t::s_nStates + dstIdx] = base_t::m_state[w * base_t::s_nStates + srcIdx];
        }
    }

    /**
     * @brief Replicate the first parallel state to all other states.
     */
    void replicateFirstState() {
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

    /**
     * @brief Orchestrate the de-phasing of parallel states using the provided jump policy.
     */
    template <typename Policy>
    void completeInitialization(const Policy& policy) {
        if (policy.hasCommon()) {
            policy.applyCommon(*this);
        }

        if constexpr (s_nStates > 1) {
            if (policy.hasSequential()) {
                policy.applySequential(*this);
            } else {
#if (RANDGEN_TESTING!=1)
                THROW("Having multiple states and no sequential jump provided does not make sense");
#endif
                replicateFirstState();
            }
        }
    }

public:
    RandGen() {}

    /**
     * @brief Template constructor for seed-based initialization.
     * Supports both Matrix and Polynomial jump arguments via parameter packs.
     */
    template <typename SeedT, typename... Args>
    RandGen(SeedT seed, Args... args) : base_t() {
        reinit(seed, args...);
    }

    /**
     * @brief Matrix-based re-initialization.
     */
    template <typename SeedT>
    void reinit(SeedT s, size_t nRepeat, const matrix_t* common, const matrix_t* sequential) {
        base_t::reinitMainState(static_cast<output_word_t>(s));
        completeInitialization(MatrixJumpPolicy{nRepeat, common, sequential});
        base_t::reinitPointers();
    }

    /**
     * @brief Polynomial-based re-initialization.
     */
    template <typename SeedT>
    void reinit(SeedT s, const poly_t* common, const poly_t* sequential) {
        base_t::reinitMainState(static_cast<output_word_t>(s));
        completeInitialization(PolyJumpPolicy{common, sequential});
        base_t::reinitPointers();
    }

    /**
     * @brief Array-based re-initialization (uint32/uint64 seeds).
     */
    template <typename SeedT>
    void reinit(const SeedT* seeds, uint32_t nSeeds, size_t nRepeat, const matrix_t* common, const matrix_t* sequential) {
        base_t::reinitMainState(seeds, nSeeds);
        completeInitialization(MatrixJumpPolicy{nRepeat, common, sequential});
        base_t::reinitPointers();
    }

    template <typename SeedT>
    void reinit(const SeedT* seeds, uint32_t nSeeds, const poly_t* common, const poly_t* sequential) {
        base_t::reinitMainState(seeds, nSeeds);
        completeInitialization(PolyJumpPolicy{common, sequential});
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

// V-Family: Inter-state vectorized generators.
template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::s_isa
         >
struct VMT19937 : details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<32>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<32>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

template < ISA Isa = SIMD_ISA
         , bool QryBlk16 = false
         >
struct XMT19937 : details::RandGen<details::MT19937Base<IsaTraits<Isa>::s_hwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<IsaTraits<Isa>::s_hwBitLen, Isa, true, QryBlk16, details::MT19937Params<32>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

// V-Family: Inter-state vectorized 64-bit MT.
template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::s_isa
         >
struct VMT19937_64 : details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<VRegBitLen, Isa, false, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

// X-Family: Intra-state vectorized 64-bit MT.
template < ISA Isa = SIMD_ISA
         , bool QryBlk16 = false
         >
struct XMT19937_64 : details::RandGen<details::MT19937Base<(IsaTraits<Isa>::s_hwBitLen >= 64 ? IsaTraits<Isa>::s_hwBitLen : 64), Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>
{
    using base_t = details::RandGen<details::MT19937Base<(IsaTraits<Isa>::s_hwBitLen >= 64 ? IsaTraits<Isa>::s_hwBitLen : 64), Isa, true, QryBlk16, details::MT19937Params<64>>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

// V-Family: Inter-state vectorized SFMT.
template < size_t VRegBitLen = SIMD_N_BITS
         , bool QryBlk16 = false
         , ISA Isa = details::BestIsa<VRegBitLen>::s_isa
         >
struct VSFMT19937 : details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>
{
    using base_t = details::RandGen<details::SFMT19937Base<VRegBitLen, Isa>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

// X-Family: Intra-state vectorized SFMT.
template < ISA Isa = SIMD_ISA
         , bool QryBlk16 = false
         >
struct XSFMT19937 : details::RandGen<details::SFMT19937Base<128, Isa>, QryBlk16>
{
    using base_t = details::RandGen<details::SFMT19937Base<128, Isa>, QryBlk16>;
    using base_t::RandGen; // reuse constructors
};

} // namespace xvmt
