# MS-MT19937

# A SIMD friendly pseudo random number generator, combining multiple Mersenne Twister 19937 generators

This repository implements a new SIMD-friendly random number generator, which maintains the same statistical properties and period as the MT19937, the well known generator proposed by Nishimura and Matsumoto.
It combines the random streams of multiple MT19937 instances, each with state vectors de-phased via jump-ahead transformations, then polls each instance in a round-robin fashion.
By evolving their vector states simultaneously, the new generator achieves perfect vectorization, fully leveraging on SIMD hardware capabilities.
Comprehensive test results demonstrate that the throughput of the new generator scales approximately linearly with the width of the SIMD registers used.

A paper describing in detail the implementation is available on [arXiv]{http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c}

## Requirements
The library is written in C++17, although with a bit of work it could be downgraded to C++03.
It is implemented with SIMD instructions available on modern X86-64 processors.
It is possible to compile it for regsiter if 128, 256 or 512 bits length.

## Library organization
The library is header only. The only source code files you need are the ones in the `include` subfolder.
In will need also the jump matrices. These are available in 7z format in the `dat` subfolder.
Only one of the matrices is needed, depending on the length of the SIMD register avialble instruction set available.
In particular, you will need one of the following files:
- 128-bits: F19935.bits
- 256-bits: F19934.bits
- 512-bits: F19933.bits
You can either uncompress the 7z files manually, or install 7za and use the `make` command.

Tested with g++11 and VS2017.

## Usage example (128 bits)
The following examples uses a 128-bits architecure is available in the file `src/demo.cpp`

```c++
    const static size_t VecLen = 128;

    // With 128 bits we use M=128/32=4 state vectors
    // We use a jump ahead matrix of 2^1993X values, where X=log2(M)
    // With M=4 we have X=5, i.e. we use the matrix for jump ahead of 2^19935 values
    // When initializing with this matrix:
    // - the 1st state is the same as the original MT19937
    // - the 2nd state is shifted forward by 1*2^19935 values
    // - the 3rd state is shifted forward by 2*2^19935 values
    // - the 4th state is shifted forward by 3*2^19935 values
    const MT19937Matrix *jumpMatrix = new MT19937Matrix("./dat/F19935.bits");

    // Optionally, we could pass an additional jump matrix affecting all states
    // This could be useful to create multiple indepdnent streams for parallelization
    // In this example we do not pass any matrix.
    const MT19937Matrix *commonJumpMatrix = nullptr;

    // This is the initialization seed. Refer to the original MT19937 documentation.
    const uint32_t seedlength = 4;
    const uint32_t seedinit[seedlength] = { 0x123, 0x234, 0x345, 0x456 };
    
    // Create the generator setting QueryMode as Block16
    MSMT19937<VecLen, QM_Block16> mt(seedinit, seedlength, commonJumpMatrix, jumpMatrix);

    // The jump matrix is no longer needed and can be released here.
    delete jumpMatrix;

    // create storage vector aligned with cache lines, where we will store results
    uint32_t* buffer = myAlignedNew<uint32_t, 64>(16);

    // Query 10 times the generator in blocks of 16 numbers
    // We could also use mt.genrand_uint32(), which queries one number at a time, but it is slower.
    for (size_t i = 0; i < 10; ++i) {
        mt.genrand_uint32_blk16(buffer);
        for (size_t j = 0; j < 16; ++j)
            std::cout << buffer[j] << ", ";
    }
    std::cout << "\n";

    // release memory
    myAlignedDelete(buffer);
```
If we were using 256 bits or 512 bits, we would only need to change the value assigned to `VecLen` and the name of the jump matrix.

## Build tests
You must first de-compress the 7z files in the dat directory.

```
** compile with SSE2 code **
make NBITS=128

** compile with AVX code **
make NBITS=256

** compile with AVX512 code **
make NBITS=512
```

## Build TestU01 tests

## Performance Stats

Performance stats obtained with a 12th Gen Intel I9-12900H CPU.
The table below shows the time in seconds to generate 3,000,000 of uniform discrete 32bit random numbers in [0,2^32-1).

```
--------------------------------------------------
SIMD Set       |  NO SIMD |   SSE4   |    AVX2   |
--------------------------------------------------
Vector Length  |     1    |     2    |      4    |
Time (seconds) |   8.4s   |   4.7s   |    4.07s  |
--------------------------------------------------
```

## TestU01

./configure --prefix=/workspace/repos/testu01/install/ --disable-shared
make -j4
make install -j4

