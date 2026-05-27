Model: SFMT19937
State Size: 19,937 bits

This directory contains jump-ahead data for the SIMD-oriented Fast Mersenne Twister.

1. Polynomial Jump Masks (New, Memory Efficient)
-----------------------------------------------
Format: Raw binary coefficients of the jump polynomial g(x) = x^J mod P(x).
Layout: Coefficients g_0, g_1, ... packed into 32-bit words (LSB first).
File Size: 4,096 bytes (fixed size for 32,768 bits max capacity).
Naming: `J<N>.sfmt.bits` represents a jump of 2^N blocks.
Usage: These masks are used by the new `PolynomialJumpApplier` and are the default for VSFMT19937.

2. Matrix Jump Files (Legacy)
-----------------------------
Format: Raw binary transition matrix F^J.
Dimensions: 19,937 x 19,937 bits.
Total File Size: ~50 MB.
Naming: `F<N>.sfmt.bits` represents a jump of 2^N blocks.
Note: Legacy matrix files are no longer stored in this branch. They are available on the `jump-matrices` branch.

3. Characteristic Polynomial
----------------------------
File: `characteristic.sfmt.hex`
Format: Hexadecimal string representing the generator's characteristic polynomial P(x).
Usage: Required by `jump_poly_generator.exe` to calculate new jump masks.
