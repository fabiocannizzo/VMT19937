#include "jump_matrix.h"

#include "../include/MSMT19937.h"

/*
    VecLen:
        must be one of {32, 128, 256, 512 }
        determines the set of SIMD instructions used and the dimension of the state vector
          32 => no SIMD, state size = 640
          128 => SSE2, state size = 640*4
          256 => SSE2, state size = 640*8
          512 => SSE2, state size = 640*16
    QueryBlkSize:
        random numbers are awlays queried in blocks of QueryBlkSize
        must be one of {1, 16, stateSize}
        16 is usually a good choice
*/
void demo128()
{

    const static size_t VecLen = 128;
    const static size_t s_M = VecLen / 32;
    std::cout << "example of how to use the generator with M=4 and SSE2 instructions\n";

    // With 128 bits we use M=128/32=4 state vectors
    // We use a jump ahead matrix of 2^1993X values, where X=log2(M), i.e. 2
    // This is the jump ahead matrix 2^19935 values
    // When initializing with this matrix:
    // - the 1st state is the same as the original MT19937
    // - the 2nd state is shifted forward by 1*2^19935 values
    // - the 3rd state is shifted forward by 2*2^19935 values
    // - the 4th state is shifted forward by 3*2^19935 values
    MT19937Matrix jumpMatrix("./dat/F19935.bits");

    // Optionally, we could pass an additional jump matrix affecting all states
    // This could be useful to create multiple indepdnent streams for parallelization
    MT19937Matrix *commonJumpMatrix = nullptr;

    const uint32_t seedlength = 4;
    const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };
    MSMT19937<128> mt(seedinit, seedlength, commonJumpMatrix, &jumpMatrix);

    // create storage vector aligned with cache lines, where we will store results
    uint32_t* buffer = myAlignedNew<uint32_t, 64>(16);

    for (size_t i = 0; i < 10; ++i) {
        mt.genrand_uint32_blk16(buffer);
        for (size_t j = 0; j < 16; ++j)
            std::cout << buffer[j] << ", ";
    }
    std::cout << "\n";

    myAlignedDelete(buffer);
}

int main()
{
    demo128();
    return 0;
}
