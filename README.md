# MS-MT19937

# A SIMD friendly pseudo random number generator, combining multiple Mersenne Twister 19937 generators

This repository implements a new SIMD-friendly random number generator, which maintains the same statistical properties and period as the MT19937, the well known generator proposed by Nishimura and Matsumoto.
It combines the random streams of multiple MT19937 instances, each with state vectors de-phased via jump-ahead transformations, then polls each instance in a round-robin fashion.
By evolving their vector states simultaneously, the new generator achieves perfect vectorization, fully leveraging on SIMD hardware capabilities.
Comprehensive test results demonstrate that the throughput of the new generator scales approximately linearly with the width of the SIMD registers used.

A paper describing in detail the implementation is available on [arXiv]{http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c}

## Requirements
The library is written in C++17. It uses extenseively Intel SIMD instructions, available on modern X86-64 processors.
It can be compiled for hardware with register length of 128, 256 or 512 bits.

## Usage
The library is header only. You only need to include the header file `include/MSMT19937.h`.
You will need also the _jump-ahead_ matrices. These are available in _7z_ format in the `dat` subfolder.
You can either uncompress the _7z_ files manually, or install `7za` and use the `make` command.
For example, tp extract the matrix `F19935.bits`, on Linux or Cygwin you can type the command `make dat/F19935.bits`.
Only one of the matrices is needed. Which one depends on the length of the SIMD registers available:
- 128-bits: `F19935.bits`
- 256-bits: `F19934.bits`
- 512-bits: `F19933.bits`

If you want to define multiple independent generator, for example to work with parallel Monte Carlo, you can use the _jump-ahead_ functionality
described later. In that case you may need to extract further jump matrices, e.g. `F00100.bits`.

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
The table below shows the time in seconds to generate 5 billions of uniform discrete 32-bit random numbers in the range $[0,2^{32}-1]$.
`NBITS` and `QUERY` are the template parameters of the generator. Peformance is compared against the original MT19937 generator and the SFMT19937 variation.
```
----------------------------------------------------------------
| GENERATOR   | NBITS | QUERY | TARGET | CPU-1 | CPU-2 | CPU-3 |
----------------------------------------------------------------
| MT19937     | n.a.  | 1    | SSE2    | 31.56 | 20.07 | 16.90 | 
| SFMT19937   | n.a.  | 1    | SSE2    | 21.67 | 6.99  | 9.97  |
| MS-MT19937  | 32    | 1    | SSE2    | 20.83 | 11.10 | 13.54 |
| MS-MT19937  | 128   | 1    | SSE2    | 13.28 | 6.19  | 7.14  |
| MS-MT19937  | 128   | 16   | SSE2    | 7.77  | 3.59  | 4.19  |
| MS-MT19937  | 128   | 2496 | SSE2    | 7.42  | 3.37  | 4.59  |
| MS-MT19937  | 256   | 1    | AVX     | n.a.  |5.43   | 6.42  |
| MS-MT19937  | 256   | 16   | AVX     | n.a.  |2.15   | 2.15  |
| MS-MT19937  | 256   | 4992 | AVX     | n.a.  |2.10   | 2.06  |
| MS-MT19937  | 512   | 1    | AVX512  | n.a.  | n.a.  | 5.66  |
| MS-MT19937  | 512   | 16   | AVX512  | n.a.  | n.a.  | 1.45  |
| MS-MT19937  | 512   | 9984 | AVX512  | n.a.  | n.a.  | 1.14  |
----------------------------------------------------------------
```
Performance stats obtained with the following CPUs:
- CPU-1: Intel(R) Celeron(R) J4125, cache 4Mb, frequency 2.0 GHz, burst frequency 2.7 GHz, SIMD support for SSE4.2. This is a low end CPU.
- CPU-2: Intel® Core™ i9-12900H, cache 24Mb cache, base frequency 3.8GHz, turbo frequency 5.0 GHz, SIMD support for AVX2. This is a high performance modern laptop CPU.
- CPU-3: Intel® Xeon® Gold 6234, cache 24.75Mb, base frequency 3.3GHz, turbo frequency 4.0 GHz, SIMD support for AVX512. This is a high performance modern desktop CPU.

## TestU01

./configure --prefix=/workspace/repos/testu01/install/ --disable-shared
make -j4
make install -j4

