# MT19937-SIMD

# A SIMD based implementation of the Mersenne Twister 19937

This generator produces a numerical sequence identical to the original MT19937 generator from Nishimura and Matsumoto:
[mt19937ar.c]{http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/mt19937ar.c}

It is twice as fast on CPUs supporting SIMD instructions.

It is a header only library

Tested with g++11 and VS2017.

## Example
```c++
   MT19937<4> mt(5489);

   for (size_t i = 0; i < 10; ++i )
       std::cout << mt.genrand_uint32() << "\n";
```

## Build
```
** compile with SSE2 code **
g++ -msse4.1 -O3 demo.cpp

** compile with AVX code **
g++ -mavx2 -O3 demo.cpp
```

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
