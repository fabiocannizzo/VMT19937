#define _CRT_SECURE_NO_WARNINGS

#include "MT19937.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>

int main()
{
    int32_t i;
    uint32_t init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    MT19937 mt(init, length);

    std::cout << "1000 outputs of genrand_uint32()\n";

    //char c;
    //std::cin >> c;

    for (i = 0; i < 1000; i++)
        std::cout << std::setw(10) << mt.genrand_uint32() << (i % 5 == 4 ? "\n" : " ");

    auto start = std::chrono::system_clock::now();
    for (size_t i = 0; i < 1000000000; ++i)
        mt.genrand_uint32();
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
        << "elapsed time: " << elapsed_seconds.count() << "s"
        << std::endl;

    return 0;
}