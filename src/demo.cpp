#include "MSMT19937.h"

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
    Use the following jump matrices, depending on the chosen number of bits:
        32  => nullptr
        128 => ./dat/F19935.bits
        256 => ./dat/F19934.bits
        512 => ./dat/F19933.bits
        
*/
void demo128()
{
    const static size_t VecLen = 128;
    std::cout << "example of how to use the generator with M=4 and SSE2 instructions\n";

    // With 128 bits we use M=128/32=4 state vectors
    // We use a jump ahead matrix of 2^1993X values, where X=log2(M)
    // With M=4 we have X=5, i.e. we use the matrix for jump ahead of 2^19935 values
    // When initializing with this matrix:
    // - the 1st state is the same as the original MT19937
    // - the 2nd state is shifted forward by 1*2^19935 values
    // - the 3rd state is shifted forward by 2*2^19935 values
    // - the 4th state is shifted forward by 3*2^19935 values
    MT19937Matrix *jumpMatrix = new MT19937Matrix(std::string("./dat/F19935.bits"));

    // Optionally, we could pass an additional jump matrix affecting all states
    // This could be useful to create multiple independent streams for parallelization
    // In this example we do not pass any matrix.
    MT19937Matrix *commonJumpMatrix = nullptr;

    // This is the initialization seed. Refer to the original MT19937 documentation.
    const uint32_t seedlength = 4;
    const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };

    // Create the generator
    MSMT19937<VecLen, QM_Block16> mt(seedinit, seedlength, commonJumpMatrix, jumpMatrix);

    // The jump matrix is no longer needed and can be released here.
    delete jumpMatrix;

    // Create storage vector aligned with cache lines, where we will store results
    uint32_t* buffer = myAlignedNew<uint32_t, 64>(16);

    // Query 10 times the generator in blocks of 16 numbers
    // We could also use mt.genrand_uint32(), which queries one number at a time, but it is slower.
    for (size_t i = 0; i < 10; ++i) {
        mt.genrand_uint32_blk16(buffer);
        for (size_t j = 0; j < 16; ++j)
            std::cout << buffer[j] << ", ";
    }
    std::cout << "\n";

    // release buffer
    myAlignedDelete(buffer);
}

int main()
{
    demo128();
    return 0;
}
