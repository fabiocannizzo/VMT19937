# MT19937-SIMD
A SIMD based implementation of the Mersenne Twister 19937

It is a header only library

Tested with g++12

Example:

```c++
   MT19937<4> mt(5489);

   for (size_t i = 0; i < 10; ++i )
       std::cout << mt.genrand_uint32() << "\n";
```

Performance stats obtained on XEON.
The table below shows the time in seconds to generate 3,000,000 of uniform discrete 32bit random numbers in [0,2^32-1).

SIMD Set       |  NO SIMD |   SSE4   |    AVX2   |
Vector Length  |     1    |     2    |      4    |
Time (seconds) |          |          |           |

