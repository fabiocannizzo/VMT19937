# VMT19937 - Vectorized Mersenne Twister with Polynomial Jump-Ahead

VMT19937 is a high-performance, SIMD-optimized implementation of the MT19937 and SFMT19937 PRNGs. It leverages jump-ahead transformations to achieve near-perfect vectorization (128/256/512-bit) by running multiple states in parallel.

## Features
- **SIMD Optimized**: Supports SSE4.2, AVX2, AVX-512, and ARM NEON.
- **Polynomial Jump-Ahead**: Efficient state advancement using compact (~4 KB) polynomial bitmasks instead of large (~50 MB) transition matrices.
- **Header-Only**: The core generator logic is implemented as a C++20 header-only library.
- **Vectorized Exponentiation**: SIMD-accelerated polynomial squaring and tree-reduction parallelization for fast jump mask generation.

## Jump Methodology
Linear Feedback Shift Registers (LFSRs) and related generators like MT19937 can be transitioned by $J$ steps by applying a jump-ahead matrix $F^J$ to the state vector $v$: $v_{t+J} = F^J v_t$. Alternatively, this can be viewed as multiplication by a polynomial $x^J$ in the ring of polynomials over GF(2) modulo the characteristic polynomial $P(x)$ of the generator: $v_{t+J} = (x^J \pmod{P(x)}) \cdot v_t$.

The polynomial method is significantly more memory-efficient as it only requires storing the $O(N)$ coefficients of the jump polynomial rather than the $O(N^2)$ entries of a transition matrix.

### References
- [1] Haramoto, H., Matsumoto, M., L'Ecuyer, P. (2008). "A Fast Jump Ahead Algorithm for Linear Recurrences in a Polynomial Space". SETA 2008.
- [2] Haramoto, H., Matsumoto, M., Nishimura, T., Panneton, F., L'Ecuyer, P. (2008). "Efficient Jump Ahead for F2-Linear Random Number Generators". INFORMS Journal on Computing.
- [3] Berlekamp, E. R. (1967). "Nonbinary BCH decoding".
- [4] Massey, J. L. (1969). "Shift-register synthesis and BCH decoding".

## Utilities

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
Calculates the jump transition matrix $F^J$ for a given generator.
```powershell
.\bin-128-cl\compute_jump_matrix.exe -g=<mt32|mt64|sfmt> [-t=<exponents>] [-p=<output_dir>]
```

## Getting Started
See `BUILD.txt` for detailed compilation and usage instructions.

### Quick Build (MSVC)
```powershell
cmd /c '"D:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" && make CXX=cl'
```

### Running Tests
```powershell
.\bin-128-sse42-cl\test.exe
```
The test suite now verifies both the legacy matrix (`F*.<gen>.bits`) and the new polynomial jump (`J*.<gen>.bits`) methods for equivalence.
