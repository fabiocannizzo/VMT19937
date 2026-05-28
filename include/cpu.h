#pragma once
#include <string>
#include <map>

// -------------------------------------------------------------
// CPUID-based info (x86/x64 only)
// -------------------------------------------------------------
struct CpuInfo
{
    std::string m_vendor;
    std::string m_brand;
    std::string m_simd;
    int m_family = 0;
    int m_model = 0;
    int m_stepping = 0;
    int m_l1d_kb = 0;
    int m_l2_kb = 0;
    int m_l3_kb = 0;
    double m_mhz = 0.0;

    std::map<std::string, size_t> m_cache;

    void print() const;
};

CpuInfo detectCpuInfo();

void setCpuAffinity(int cpu);
void setPriorityHigh();

