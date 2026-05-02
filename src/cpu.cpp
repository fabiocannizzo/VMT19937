#include "../include/cpu.h"

#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <cstring>
#include <chrono>
#include <fstream>

#if defined(_WIN32) || defined(__CYGWIN__)
    #include <windows.h>
    #include <intrin.h>
#else
    #include <sys/resource.h>
    #if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
        #include <cpuid.h>
    #endif
    #include <sched.h>
    #include <unistd.h>
    #include <fstream>
    #include <sstream>
#endif


// simple helper for CPUID
static inline void cpuid(int regs[4], int eax, int ecx = 0)
{
#if defined(_WIN32)
    __cpuidex(regs, eax, ecx);
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    __cpuid_count(eax, ecx, regs[0], regs[1], regs[2], regs[3]);
#else
    // Not an x86 architecture
    regs[0] = regs[1] = regs[2] = regs[3] = 0;
#endif
}

CpuInfo detectCpuInfo()
{
    CpuInfo info;

#if !defined(__x86_64__) && !defined(_M_X64) && !defined(__i386__) && !defined(_M_IX86) && !defined(_WIN32)
    // ARM/Generic Linux detection
    info.vendor = "ARM";
    info.simd = "NEON";

    std::ifstream f("/proc/cpuinfo");
    std::string line;
    while (std::getline(f, line))
    {
        if (line.find("Model") != std::string::npos || line.find("Hardware") != std::string::npos)
        {
            size_t pos = line.find(":");
            if (pos != std::string::npos) info.brand = line.substr(pos + 2);
        }
        if (line.find("cpu MHz") != std::string::npos || line.find("BogoMIPS") != std::string::npos)
        {
            double val;
            if (sscanf(line.c_str(), "%*s %*s : %lf", &val) == 1) info.mhz = val;
        }
    }
    return info;
#else
    int regs[4] = {0};

    // Vendor
    cpuid(regs, 0);
    char vendor[13];
    memcpy(vendor + 0, &regs[1], 4);
    memcpy(vendor + 4, &regs[3], 4);
    memcpy(vendor + 8, &regs[2], 4);
    vendor[12] = '\0';
    info.vendor = vendor;

    // Brand string
    char brand[49];
    memset(brand, 0, sizeof(brand));
    cpuid(regs, 0x80000000);
    unsigned int maxExtended = regs[0];
    if (maxExtended >= 0x80000004)
    {
        int* p = (int*)brand;
        for (int i = 0; i < 3; ++i)
        {
            cpuid(regs, 0x80000002 + i);
            memcpy(p + i * 4, regs, sizeof(regs));
        }
        info.brand = brand;
    }

    // Family, model, stepping
    cpuid(regs, 1);
    info.stepping = regs[0] & 0xF;
    info.model = (regs[0] >> 4) & 0xF;
    info.family = (regs[0] >> 8) & 0xF;

    // SIMD detection
    bool sse2 = regs[3] & (1 << 26);
    bool sse3 = regs[2] & (1 << 0);
    bool ssse3 = regs[2] & (1 << 9);
    bool sse41 = regs[2] & (1 << 19);
    bool sse42 = regs[2] & (1 << 20);
    bool avx = regs[2] & (1 << 28);

    cpuid(regs, 7, 0);
    bool avx2 = regs[1] & (1 << 5);
    bool avx512f = regs[1] & (1 << 16);

    if (avx512f)
        info.simd = "AVX-512";
    else if (avx2)
        info.simd = "AVX2";
    else if (avx)
        info.simd = "AVX";
    else if (sse42)
        info.simd = "SSE4.2";
    else if (sse41)
        info.simd = "SSE4.1";
    else if (ssse3)
        info.simd = "SSSE3";
    else if (sse3)
        info.simd = "SSE3";
    else if (sse2)
        info.simd = "SSE2";
    else
        info.simd = "None";

#ifdef _WIN32
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    DWORD bufferSize = 0;
    DWORD mhz = 0;
    bufferSize = sizeof(DWORD);
    if (RegGetValueA(HKEY_LOCAL_MACHINE,
                     "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
                     "~MHz",
                     RRF_RT_DWORD, nullptr, &mhz, &bufferSize) == ERROR_SUCCESS)
        info.mhz = (double)mhz;
#else
    std::ifstream f("/proc/cpuinfo");
    std::string line;
    while (std::getline(f, line))
    {
        if (line.find("cache size") != std::string::npos)
        {
            int size;
            char unit[8];
            sscanf(line.c_str(), "cache size\t: %d %s", &size, unit);
            info.l2_kb = size; // approximate
        }
        if (line.find("cpu MHz") != std::string::npos)
        {
            double val;
            sscanf(line.c_str(), "cpu MHz\t\t: %lf", &val);
            info.mhz = val;
        }
    }
#endif

    for (int i = 0; ; ++i) {
        cpuid(regs, 4, i);;
        unsigned cacheType = regs[0] & 0x1F;
        if (cacheType == 0) break; // No more caches

        unsigned level = (regs[0] >> 5) & 0x7;
        unsigned sets = regs[2] + 1;
        unsigned ways = ((regs[1] >> 22) & 0x3FF) + 1;
        unsigned partitions = ((regs[1] >> 12) & 0x3FF) + 1;
        unsigned lineSize = (regs[1] & 0xFFF) + 1;
        unsigned cacheSize = sets * ways * partitions * lineSize;

        // Determine cache type
        const char* typeStr = nullptr;
        switch (cacheType)
        {
            case 1: typeStr = "d"; break; // Data cache
            case 2: typeStr = "i"; break; // Instruction cache
            case 3: typeStr = "u"; break; // Unified cache
            default: typeStr = "?"; break;
        }

        std::ostringstream key;
        key << "L" << level << typeStr;
        info.cache[key.str()] = cacheSize;
    }
#endif

    return info;
}

// -------------------------------------------------------------
// Affinity and priority (as before)
// -------------------------------------------------------------
void setCpuAffinity(int cpu)
{
#if defined(_WIN32) || defined(__CYGWIN__)
    DWORD_PTR mask = 1ULL << cpu;
    auto ret = SetProcessAffinityMask(GetCurrentProcess(), mask);
    if (ret == 0) {
        std::cout << "WARNING: failure setting affinity\n";
    }
#elif defined(__linux__)
    cpu_set_t cpuset;
    std::memset(&cpuset, 0, sizeof(cpuset));
    ((unsigned long*)&cpuset)[0] |= 1UL << cpu;
    if (sched_setaffinity(0, sizeof(cpuset), &cpuset) != 0)
        std::cout << "WARNING: failure setting affinity\n";
#else
#error Unsupported platform
#endif
    std::cout << "Set CPU affinity to CPU " << cpu << "\n";
}

void setPriorityHigh()
{
#if defined(_WIN32) || defined(__CYGWIN__)
    SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
#else
    setpriority(PRIO_PROCESS, 0, -10);
#endif
    std::cout << "Set high process priority\n";
}

// -------------------------------------------------------------
void CpuInfo::print() const
{
    std::cout << "=== CPU Information ===\n";
    std::cout << "Vendor : " << vendor << "\n";
    std::cout << "Brand  : " << brand << "\n";
    std::cout << "Family : " << family << ", Model: " << model
              << ", Stepping: " << stepping << "\n";
    std::cout << "SIMD   : " << simd << "\n";
    std::cout << "Clock  : " << mhz << " MHz\n";
    for (auto&[k,s] : cache)
        std::cout << k << " Cache : Level " << k << ", size: " << s / 1024 << " KB\n";
    std::cout << "=======================\n";
}
