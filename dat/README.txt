Jump Data Organization
----------------------

This directory contains data for state advancement (jump-ahead):

1. dat/poly/   - Polynomial jump-ahead bitmasks (J*.bits) and characteristic polynomials (.hex).
                 These are compact (~4 KB) and part of the main branch.
2. dat/matrix/ - Legacy jump transition matrices (F*.bits).
                 These are not included in this branch to keep the repo size small.
                 To retrieve them, run: git checkout jump-matrix -- dat/matrix/
                 Then decompress any .7z files using 7-Zip.

Generating Data
---------------
Jump files can be regenerated from scratch using the jump-ahead utilities.
See BUILD.txt in the repository root for detailed instructions.
