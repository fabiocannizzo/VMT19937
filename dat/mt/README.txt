Model: MT19937
Dimensions: 19,937 rows x 19,937 columns
Storage Format: Raw binary.
Layout: Row-major order. Each row is stored as a contiguous sequence of bits packed into bytes (8 bits per byte, MSB of first byte is first bit of row).
Row Size: 2,493 bytes (19,937 bits / 8, rounded up).
Total File Size: 49,702,941 bytes.
File Naming Convention: `F<N>.bits` contains the transition matrix for a jump of size 2^N.

Note: The actual binary matrix files (.bits) are not stored in this branch to keep the repository size manageable. They must be obtained by checking out the `jumpfiles` branch.
