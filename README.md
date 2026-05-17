# VMT19937

A high-performance C++20 library for SIMD-vectorized Mersenne Twister 19937 PRNGs, supporting x86-64 (SSE4.2, AVX2, AVX-512) and ARM64 (NEON, SVE).

## Introduction

VMT19937 provides two families of SIMD-optimized generators based on the Mersenne Twister 19937:

1.  **X-Family (Intra-state Vectorization):** Vectorizes the state recurrence of a single MT instance. It maintains the exact bit-for-bit mathematical identity and sequential order of the original algorithms.
    *   `XMT19937`: Vectorized 32-bit MT19937 [^1].
    *   `XMT19937_64`: Vectorized 64-bit MT19937-64 [^2].
    *   `XSFMT19937`: SIMD-oriented Fast Mersenne Twister [^3].
2.  **V-Family (Inter-state Vectorization):** Combines multiple independent MT instances, de-phased via jump-ahead transformations, and polls them in round-robin fashion. This approach achieves perfect vectorization and scales linearly with SIMD register width.
    *   `VMT19937`: Multi-state 32-bit MT19937.
    *   `VMT19937_64`: Multi-state 64-bit MT19937-64.
    *   `VSFMT19937`: Multi-state SIMD-oriented Fast Mersenne Twister.

The library is header-only and leverages modern C++20 features and architecture-specific intrinsics (including AVX-512 ternary logic and SVE/NEON bridges) to outperform proprietary vendor libraries like Intel MKL.

## Parametrization

Generators are template-parameterized to allow fine-tuning of performance and ISA compatibility:

### 1. V-Family (Inter-state Vectorization)
The V-Family uses multiple parallel states to achieve full SIMD width utilization.

```cpp
template <
    size_t VRegBitLen = SIMD_N_BITS,  // Logical SIMD width (128, 256, 512)
    bool QryBlk16 = false,            // Toggle for optimized Block-16 query mode
    ISA Isa = details::BestIsa<VRegBitLen>::isa // Target ISA (SSE42, AVX2, AVX512, NEON, SVE)
>
class VMT19937;
```

*   **`VRegBitLen`**: The logical register width. This parameter determines the **mathematical sequence** (i.e., the number of parallel states).
*   **`Isa`**: Determines the hardware **optimization level**.
*   **Portability Remark:** By fixing `VRegBitLen`, you ensure that the generator produces the exact same sequence regardless of the underlying hardware ISA. For example, `VMT19937<512>` will produce the same stream on an AVX2 machine (where it is emulated) and an AVX-512 machine (where it is native). `VRegBitLen` can be larger than the hardware register size (`HwSize`), allowing for seamless portability across different architectures. In a nutshell, the sequence is decided by `VRegBitLen`, the optimization by the `Isa`.
*   **`QryBlk16`**: When `true`, enables the optimized `genrand_uint32_blk16()` and `genrand_word_blk()` methods.

### 2. X-Family (Intra-state Vectorization)
The X-Family vectorizes the internal recurrence of a single MT state. It maintains sequential consistency with the original MT19937 algorithm.

```cpp
template <
    ISA Isa = SIMD_ISA,               // Target ISA (SSE42, AVX2, AVX512, NEON, SVE)
    bool QryBlk16 = false             // Toggle for optimized Block-16 query mode
>
class XMT19937;
```

*   **`Isa`**: Determines the hardware instruction set used to accelerate the single-state recurrence.
*   **`QryBlk16`**: Same as in the V-Family, enables zero-overhead bulk generation.
*   *Note: X-Family generators automatically use the full hardware register width for the target ISA.*

## Jump Matrices

For the **V-Family** generators to produce independent streams, the internal parallel states must be initialized with a **jump matrix**. The choice of the matrix depends on the number of states $N$, which is determined by the logical register width (`VRegBitLen`).

The number of states is calculated as $N = VRegBitLen / StateWordSize$, where **StateWordSize** is:
*   **32** for `VMT19937`
*   **64** for `VMT19937_64`
*   **128** for `VSFMT19937`

### Matrix Selection Formula
The goal is to partition the period $P = 2^{19937}-1$ into $N$ equal segments of length $L = P/N \approx 2^{19937-k}$, where $k = \log_2(N)$. To ensure maximum separation between the parallel states, you should use the matrix $F_{19937-k}$ (found in the `dat/` folder).

### Recommended Matrices

| Generator | StateWordSize | Width (`VRegBitLen`) | States ($N$) | Recommended Matrix |
| :--- | :---: | :---: | :---: | :--- |
| **VMT19937** | 32 | 128 | 4 | `dat/mt32/F19935.bits` |
| | 32 | 256 | 8 | `dat/mt32/F19934.bits` |
| | 32 | 512 | 16 | `dat/mt32/F19933.bits` |
| **VMT19937_64** | 64 | 128 | 2 | `dat/mt64/F19936.bits` |
| | 64 | 256 | 4 | `dat/mt64/F19935.bits` |
| | 64 | 512 | 8 | `dat/mt64/F19934.bits` |
| **VSFMT19937** | 128 | 256 | 2 | `dat/sfmt/F19936.bits` |
| | 128 | 512 | 4 | `dat/sfmt/F19935.bits` |

*Note: The **X-Family** (single-state) does not require a jump matrix for internal initialization, as it vectorizes the recurrence of a single state.*

## Performance Summary

Throughput measured in **Million samples per second (M/s)**. Benchmarks performed on Intel Celeron (SSE4.2), Xeon Platinum 8375C (AVX-512), Xeon E5-2686 v4 (AVX2), Neoverse-N1 (NEON), and Neoverse-V1 (SVE256).

### MT19937 Family (generates 32-bit random numbers)
| Mode | Generator | SSE4.2 | AVX2 | AVX-512 | NEON | SVE256 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Scalar** | ORIG-MT19937 | 161 | 259 | 418 | 227 | 362 |
| | STL-MT19937 | 134 | 306 | 599 | 277 | 570 |
| | MKL-MT* | 23 | 41 | 53 | n.a. | n.a. |
| | **X-MT19937** | 317 | 417 | **1341** | 421 | **672** |
| | **V-MT19937** | **341** | **448** | 729 | **431** | 652 |
| **Vectorial** | ORIG-MT19937??? | 161 | 259 | 418 | 227 | 362 |
| | STL-MT19937 | 135 | 303 | 565 | 279 | 548 |
| | MKL-MT | 711 | 1668 | **6479** | n.a. | n.a. |
| | **X-MT19937** | 711 | 1969 | 6418 | 583 | 1618 |
| | **V-MT19937** | **858** | **2392** | 5742 | **634** | **1869** |

### SFMT19937 Family (generates 32-bit random numbers)
| Mode | Generator | SSE4.2 | AVX2 | AVX-512 | NEON | SVE256 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Scalar** | ORIG-SFMT | 440 | 792 | 1180 | 524 | 752 |
| | MKL-SFMT* | 21 | 38 | 61 | n.a. | n.a. |
| | **X-SFMT19937** | **472** | 848 | 1193 | **555** | 672 |
| | **V-SFMT19937** | n.a. | **944** | **1507** | n.a. | **898** |
| **Vectorial** | ORIG-SFMT | 1444 | 2504 | 4424 | **1381** | 1692 |
| | MKL-SFMT | 1654 | 3360 | 2761 | n.a. | n.a. |
| | **X-SFMT19937** | **1700** | 2308 | 3584 | 1234 | 1327 |
| | **V-SFMT19937** | n.a. | **4515** | **10961** | n.a. | **2960** |

### MT19937-64 Family (generates 64-bit random numbers)
| Mode | Generator | SSE4.2 | AVX2 | AVX-512 | NEON | SVE256 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Scalar** | ORIG-MT-64 | 140 | 180 | 338 | 203 | 337 |
| | STL-MT-64 | 165 | 317 | 590 | 285 | **573** |
| | **X-MT19937-64** | 180 | 325 | 680 | 284 | 556 |
| | **V-MT19937-64** | **208** | **372** | **851** | **334** | 513 |
| **Vectorial** | ORIG-MT-64??? | 140 | 180 | 338 | 203 | 337 |
| | STL-MT-64 | 165 | 314 | 565 | 281 | 542 |
| | **X-MT19937-64** | 252 | 827 | 1547 | 255 | **927** |
| | **V-MT19937-64** | **378** | **1161** | **2886** | **290** | 901 |

<sup>*</sup>MKL tested in vectorial mode with block size 1.
<sup>???</sup>No native vectorial support; values for ORIG-* copied from scalar mode; STL-* performance reflects compiler auto-vectorization.

### STL Auto-Vectorization
An interesting phenomenon was observed with the C++ Standard Template Library (STL) implementation of MT19937. While the STL API does not provide a vectorial interface, modern compilers like GCC are capable of auto-vectorizing the internal tempering logic when the generator is polled in a tight loop.

As shown in the performance tables, the STL generators exhibit significantly higher throughput when the compiler is free to parallelize the operations, often matching or exceeding the performance of explicit SIMD implementations on certain architectures (e.g., AVX-512 and SVE256). To ensure a fair comparison in scalar mode, all benchmarks were updated to use a `volatile` variable trap, which prevents the compiler from using Dead-Code Elimination (DCE) or auto-vectorization, thus reflecting true sequential execution.

## Usage

### 1. Basic Scalar Usage
```cpp
#include "RandGen.h"

// XMT19937 is single-state and doesn't require jump matrices.
xvmt::XMT19937<> gen(42); // Initialize with seed 42
uint32_t val = gen.genrand_uint32();
```

### 2. High-Performance Bulk Generation (V-Family)
The V-Family requires a **jump matrix** to initialize de-phased states. Matrices are provided in the `dat/` folder.

```cpp
#include "RandGen.h"
#include <memory>

using namespace xvmt;

// 1. Load the appropriate jump matrix (e.g., for 256-bit AVX2)
auto jumpMatrix = std::make_unique<MT19937Matrix<32>>("./dat/mt32/F19934.bits");

// 2. Initialize the generator with Block-16 optimization
VMT19937<256, true> gen(1234, 0, nullptr, jumpMatrix.get());

// 3. Generate into cache-aligned memory
alignas(64) uint32_t buffer[16];
gen.genrand_uint32_blk16(buffer);
```

### 3. Multiple Independent Generators (Parallel Streams)
For parallel Monte Carlo or multi-threaded applications, use a **common jump matrix** to ensure each generator instance produces a non-overlapping stream.

```cpp
#include "RandGen.h"
#include <array>
#include <memory>

using namespace xvmt;

// Matrix to separate internal parallel states (SSE 128-bit)
auto seqJump = std::make_unique<MT19937Matrix<32>>("./dat/mt32/F19935.bits");

// Matrix to separate different generator instances (by 2^100 steps)
auto commonJump = std::make_unique<MT19937Matrix<32>>("./dat/mt32/F00100.bits");

// Create 10 independent generators
std::array<std::unique_ptr<VMT19937<128, true>>, 10> gens;
for (size_t i = 0; i < 10; ++i) {
    // Each instance 'i' is jumped forward by i * 2^100 steps
    gens[i] = std::make_unique<VMT19937<128, true>>(seed, 4, i, commonJump.get(), seqJump.get());
}
```

## Build Requirements
*   **Compiler:** GCC 13+, Clang 16+, or MSVC 2022 (with `/std:c++20`).
*   **Hardware:**
    *   x86-64 with SSE4.2, AVX2, or AVX-512.
    *   ARM64 with NEON or SVE/SVE2.
*   **Build System:** CMake or the provided `Makefile`.

```bash
# Example: Build benchmarks for AVX-512
make NBITS=512
./bin/perf --dir dat
```

## References
[^1]: 1998, M. Matsumoto, T. Nishimura, "Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator", ACM TOMACS.
[^2]: 2000, T. Nishimura, "Tables of 64-Bit Mersenne Twisters", ACM TOMACS.
[^3]: 2008, M. Saito and M. Matsumoto, "SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number Generator", Springer.
[^4]: 2007, P. L'Ecuyer and R. Simard, "TestU01: A C Library for Empirical Testing of Random Number Generators", ACM TOMS.

## License
MIT License. Includes reference code for MT19937 and SFMT for benchmarking purposes.

