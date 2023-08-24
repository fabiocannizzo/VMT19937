#pragma once
#include <string>
#include <map>

// -------------------------------------------------------------
// CPUID-based info (x86/x64 only)
// -------------------------------------------------------------
struct CpuInfo
{
    std::string vendor;
    std::string brand;
    std::string simd;
    int family = 0;
    int model = 0;
    int stepping = 0;
    int l1d_kb = 0;
    int l2_kb = 0;
    int l3_kb = 0;
    double mhz = 0.0;

    std::map<std::string, size_t> cache;

    void print() const;
};

CpuInfo detectCpuInfo();

void setCpuAffinity(int cpu);
void setPriorityHigh();

