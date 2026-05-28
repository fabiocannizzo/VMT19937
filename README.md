# VMT19937 - Vectorized Mersenne Twister with Polynomial Jump-Ahead

A high-performance C++20 library for SIMD-vectorized Mersenne Twister 19937 PRNGs, supporting x86-64 (SSE4.2, AVX2, AVX-512) and ARM64 (NEON, SVE). VMT19937 leverages jump-ahead transformations to achieve near-perfect vectorization (128/256/512-bit) by running multiple states in parallel.

## Features
- **SIMD Optimized**: Supports SSE4.2, AVX2, AVX-512, and ARM NEON/SVE.
- **Header-Only**: The core generator logic is implemented as a C++20 header-only library.

## Introduction

VMT19937 provides two families of SIMD-optimized generators based on the Mersenne Twister 19937:

1.  **V-Family (Inter-state Vectorization):** Combines multiple independent MT instances, de-phased via jump-ahead transformations, and polls them in round-robin fashion. This approach achieves perfect vectorization and scales linearly with SIMD register width.
    *   `VMT19937`: Multi-state 32-bit MT19937.
    *   `VMT19937_64`: Multi-state 64-bit MT19937-64.
    *   `VSFMT19937`: Multi-state SIMD-oriented Fast Mersenne Twister.
2.  **X-Family (Intra-state Vectorization):** Vectorizes the state recurrence of a single MT instance. It maintains the exact bit-for-bit mathematical identity and sequential order of the original algorithms.
    *   `XMT19937`: Vectorized 32-bit MT19937 [^1].
    *   `XMT19937_64`: Vectorized 64-bit MT19937-64 [^2].
    *   `XSFMT19937`: SIMD-oriented Fast Mersenne Twister [^3].

## Parametrization

Generators are template-parameterized to allow fine-tuning of performance and ISA compatibility:

### 1. V-Family (Inter-state Vectorization)
The V-Family uses multiple parallel states to achieve full SIMD width utilization.

```cpp
template <
    size_t VRegBitLen = SIMD_N_BITS,  // Logical SIMD width (128, 256, 512)
    bool QryBlk16 = false,            // Toggle for optimized Block-16 query mode
    ISA Isa = details::BestIsa<VRegBitLen>::isa // Target ISA
>
class VMT19937;
```

*   **`VRegBitLen`**: The logical register width. This determines the **mathematical sequence** (i.e., the number of parallel states).
*   **`Isa`**: Determines the hardware **optimization level**.
*   **Portability Remark:** By fixing `VRegBitLen`, you ensure that the generator produces the exact same sequence regardless of the underlying hardware ISA. For example, `VMT19937<512>` will produce the same stream on an AVX2 machine (where it is emulated) and an AVX-512 machine (where it is native).
*   **`QryBlk16`**: When `true`, enables the optimized `genrand_uint32_blk16()` and `genrand_word_blk()` methods.

### 2. X-Family (Intra-state Vectorization)
The X-Family vectorizes the internal recurrence of a single MT state.

```cpp
template <
    ISA Isa = SIMD_ISA,               // Target ISA
    bool QryBlk16 = false             // Toggle for optimized Block-16 query mode
>
class XMT19937;
```

## State Advancement (Jump Methodology)

For the **V-Family** generators to produce independent streams, the internal parallel states must be initialized using jump-ahead transformations. VMT19937 supports two methods for this:

1. **Polynomial Jump (Default)**: Linear Feedback Shift Registers (LFSRs) and related generators can be transitioned by $J$ steps by multiplying the state by a polynomial $x^J$ in the ring of polynomials over GF(2) modulo the characteristic polynomial $P(x)$. This requires storing only $O(N)$ coefficients (~4 KB).
2. **Matrix Jump (Used for Testing)**: Applying a transition matrix $F^J$ to the state vector $v$: $v_{t+J} = F^J v_t$. This requires storing $O(N^2)$ entries (~50 MB for SFMT). Legacy matrices are stored in a separate `jump-matrix` branch.

### Jump Selection

The number of parallel states is calculated as $N = VRegBitLen / StateWordSize$, where **StateWordSize** is 32 for `VMT19937`, 64 for `VMT19937_64`, and 128 for `VSFMT19937`.

To ensure maximum separation between the $N$ parallel states, partition the period $P = 2^{19937}-1$ into segments of length $P/N \approx 2^{19937-k}$, where $k = \log_2(N)$. Use the corresponding polynomial bitmask $J_{19937-k}$.

| Generator | Width (`VRegBitLen`) | States ($N$) | Recommended Polynomial | Legacy Matrix |
| :--- | :---: | :---: | :--- | :--- |
| **VMT19937** | 128 | 4 | `dat/poly/mt32/J19935.mt32.bits` | `dat/matrix/mt32/F19935.bits` |
| | 256 | 8 | `dat/poly/mt32/J19934.mt32.bits` | `dat/matrix/mt32/F19934.bits` |
| | 512 | 16 | `dat/poly/mt32/J19933.mt32.bits` | `dat/matrix/mt32/F19933.bits` |
| **VMT19937_64** | 128 | 2 | `dat/poly/mt64/J19936.mt64.bits` | `dat/matrix/mt64/F19936.bits` |
| | 256 | 4 | `dat/poly/mt64/J19935.mt64.bits` | `dat/matrix/mt64/F19935.bits` |
| | 512 | 8 | `dat/poly/mt64/J19934.mt64.bits` | `dat/matrix/mt64/F19934.bits` |
| **VSFMT19937** | 256 | 2 | `dat/poly/sfmt/J19936.sfmt.bits` | `dat/matrix/sfmt/F19936.bits` |
| | 512 | 4 | `dat/poly/sfmt/J19935.sfmt.bits` | `dat/matrix/sfmt/F19935.bits` |

## Utilities

Includes several tools for calculating jump data.

### `characteristic_poly_finder.exe`
Derives the characteristic polynomial $P(x)$ for a given generator using the Berlekamp-Massey algorithm.
```powershell
.\bin-128-cl\characteristic_poly_finder.exe -g=<mt32|mt64|sfmt> [-o=<output_file>.hex]
```

### `jump_poly_generator.exe`
Calculates the jump polynomial $x^J \pmod{P(x)}$ for any jump step $J$.
```powershell
.\bin-128-cl\jump_poly_generator.exe -g=<mt32|mt64|sfmt> [-n=<power_of_2> | -t=<exact_steps>] -p=<char_poly_file>.hex -o=<output_mask>.bits  
```

### `compute_jump_matrix.exe` (Legacy)
Calculates the legacy jump transition matrix $F^J$ for a given generator.
```powershell
.\bin-128-cl\compute_jump_matrix.exe -g=<mt32|mt64|sfmt> [-t=<exponents>] [-p=<output_dir>]
```

## Usage

### 1. Basic Scalar Usage
```cpp
#include "RandGen.h"

// X-Family is single-state and doesn't require jump masks.
xvmt::XMT19937<> gen(42); // Initialize with seed 42
uint32_t val = gen.genrand_uint32();
```

### 2. High-Performance Bulk Generation (Polynomial Jump)
The V-Family requires a **jump mask** to initialize de-phased states.

```cpp
#include "RandGen.h"
#include "polynomial_jump.h"
#include <fstream>
#include <memory>

using namespace xvmt;

// 1. Load the appropriate polynomial jump mask (e.g., for 256-bit AVX2)
details::Polynomial<> jumpMask;
jumpMask.fromBinFile("./dat/poly/mt32/J19934.mt32.bits");

// 2. Initialize the generator
VMT19937<256> gen(1234, nullptr, &jumpMask);

// 3. Generate numbers
uint32_t val = gen.genrand_uint32();
```

### 3. Multiple Independent Generators (Parallel Streams)
For parallel Monte Carlo or multi-threaded applications, use a **common jump mask** to ensure each generator instance produces a non-overlapping stream.

```cpp
#include "RandGen.h"
#include "polynomial_jump.h"
#include <array>
#include <fstream>
#include <memory>

using namespace xvmt;

// Mask to separate internal parallel states (SSE 128-bit)
details::Polynomial<> seqJump;
seqJump.fromBinFile("./dat/poly/mt32/J19935.mt32.bits");

// Mask to separate different generator instances (by 2^100 steps)
details::Polynomial<> commonJump;
commonJump.fromBinFile("./dat/poly/mt32/J00100.mt32.bits");

// Create 10 independent generators
std::array<std::unique_ptr<VMT19937<128>>, 10> gens;
for (size_t i = 0; i < 10; ++i) {
    // Each instance 'i' is initialized with a different common jump
    // Note: To jump different amounts using polynomials, one would typically
    // pre-calculate commonJump^i. In this simple example, they all jump by 2^100.
    gens[i] = std::make_unique<VMT19937<128>>(42, &commonJump, &seqJump);
}
```

## Performance Summary

*(Throughput measured in Million samples per second (M/s).)*

### MT19937 Family (generates 32-bit random numbers)
| Mode | Generator | SSE4.2 | AVX2 | AVX-512 | NEON | SVE256 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Scalar** | ORIG-MT19937 | 161 | 259 | 418 | 227 | 362 |
| | **X-MT19937** | 317 | 417 | **1341** | 421 | **672** |
| | **V-MT19937** | **341** | **448** | 729 | **431** | 652 |
| **Vectorial** | ORIG-MT19937 | 161 | 259 | 418 | 227 | 362 |
| | MKL-MT | 711 | 1668 | **6479** | n.a. | n.a. |
| | **X-MT19937** | 711 | 1969 | 6418 | 583 | 1618 |
| | **V-MT19937** | **858** | **2392** | 5742 | **634** | **1869** |

## Build Requirements
*   **Compiler:** GCC 13+, Clang 16+, or MSVC 2022 (with `/std:c++20`).
*   **Hardware:** x86-64 (SSE4.2, AVX2, AVX-512) or ARM64 (NEON, SVE/SVE2).
*   **Build System:** CMake or the provided `Makefile`.

See `BUILD.txt` for detailed compilation and MSVC environment instructions.

## References
[^1]: 1998, M. Matsumoto, T. Nishimura, "Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator", ACM TOMACS.
[^2]: 2000, T. Nishimura, "Tables of 64-Bit Mersenne Twisters", ACM TOMACS.
[^3]: 2008, M. Saito and M. Matsumoto, "SIMD-oriented Fast Mersenne Twister", Springer.
[^4]: Haramoto, H., et al. (2008). "A Fast Jump Ahead Algorithm for Linear Recurrences in a Polynomial Space". SETA 2008.

## License
MIT License. Includes reference code for MT19937 and SFMT for benchmarking purposes.
