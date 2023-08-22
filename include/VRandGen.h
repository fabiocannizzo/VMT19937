#pragma once

#include "VMT19937.h"
#include "VSFMT19937.h"

enum VRandGenQueryMode { QM_Scalar, QM_Block16, QM_StateSize };

namespace Details
{

template < template <size_t L> class GenBase, size_t RegisterBitLen, VRandGenQueryMode QueryMode>
class VRandGen : public GenBase<RegisterBitLen>
{
    typedef GenBase<RegisterBitLen> base_t;

public:
    static const VRandGenQueryMode s_queryMode = QueryMode;

    using base_t::base_t; // use base class constructors

    // generates a random number on [0,0xffffffff] interval
    uint32_t FORCE_INLINE genrand_uint32()
    {
        static_assert(QueryMode == QM_Scalar);
        return base_t::genrand_uint32();
    }

    // generates 16 uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_blk16(uint32_t* dst)
    {
        static_assert(QueryMode == QM_Block16, "This function can only be invoked when query mode is QM_Block16");
        base_t::genrand_uint32_blk16(dst);
    }

    // generates a block of the same size as the state vector of uniform discrete random numbers in [0,0xffffffff] interval
    // for optimal performance the vector dst should be aligned on a 64 byte boundary
    void genrand_uint32_stateBlk(uint32_t* dst)
    {
        static_assert(QueryMode == QM_StateSize, "This function can only be invoked when query mode is QM_StateSize");
        base_t::genrand_uint32_stateBlk(dst);
    }
};

} // namespace Details

template <size_t RegisterBitLen = SIMD_N_BITS, VRandGenQueryMode QueryMode = QM_Scalar>
using VMT19937 = Details::VRandGen<Details::VMT19937Base, RegisterBitLen, QueryMode>;

template <size_t RegisterBitLen = SIMD_N_BITS, VRandGenQueryMode QueryMode = QM_Scalar>
using VSFMT19937 = Details::VRandGen<Details::VSFMT19937Base, RegisterBitLen, QueryMode>;
