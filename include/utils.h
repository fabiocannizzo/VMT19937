#pragma once

#include <cstdint>

// nAlign must be a multiple of 2 and no more than 64
template <typename T, unsigned nAlign>
T* myAlignedNew(size_t n)
{
    // allocate memory
    uint8_t *p = new uint8_t[n * sizeof(T) + nAlign];

    // align pointer
    auto addr = reinterpret_cast<std::uintptr_t>(p);
    uint8_t moveFwdBy = (nAlign - (addr % nAlign));
    p += moveFwdBy;
    p[-1] = moveFwdBy;

    return (T*)p;
}

void myAlignedDelete(void *p)
{
    if (p) {
        uint8_t* p8 = (uint8_t*)p;
        p8 -= p8[-1];
        delete p8;
    }
}