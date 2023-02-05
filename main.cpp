#include "MT19937.h"

#include <iostream>
#include <iomanip>

int main()
{
    int32_t i;
    uint32_t init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    MT19937 mt(init, length);

    std::cout << "1000 outputs of genrand_uint32()\n";

    for (i = 0; i < 1000; i++)
        std::cout << std::setw(10) << mt.genrand_uint32() << (i % 5 == 4 ? "\n" : " ");

    return 0;
}