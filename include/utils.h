#pragma once

#include <cstdint>

#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <type_traits>

typedef std::map<std::string, std::string> ArgMap;

inline ArgMap parseArgs(int argc, const char** argv)
{
    ArgMap args;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.empty()) continue;
        size_t start = 0;
        if (s.size() >= 2 && s[0] == '-' && s[1] == '-') start = 2;
        else if (s.size() >= 1 && s[0] == '-') start = 1;
        else continue;

        std::string body = s.substr(start);
        size_t eq = body.find('=');
        if (eq != std::string::npos) {
            args[body.substr(0, eq)] = body.substr(eq + 1);
        }
        else {
            args[body] = "";
        }
    }
    return args;
}

template <typename T>
inline bool consumeArg(ArgMap& args, const std::string& key, bool compulsory, T& result)
{
    auto it = args.find(key);
    if (it != args.end()) {
        std::string val = it->second;
        args.erase(it);

        if constexpr (std::is_same_v<T, std::string>) {
            result = val;
        }
        else {
            if (val.empty()) {
                throw std::runtime_error("Argument " + key + " requires a value");
            }
            std::stringstream ss(val);
            if (!(ss >> result) || !ss.eof()) {
                throw std::runtime_error("Argument conversion failed for key: " + key + " with value: " + val);
            }
        }
        return true;
    }

    if (compulsory) {
        throw std::runtime_error("Compulsory argument missing: " + key);
    }

    return false;
}

inline bool consumeArg(ArgMap& args, const std::string& key)
{
    auto it = args.find(key);
    if (it != args.end()) {
        args.erase(it);
        return true;
    }
    return false;
}

inline void waitForDebugger(ArgMap& args)
{
    if (consumeArg(args, "w") || consumeArg(args, "wait")) {
        std::cout << "press a key to continue..." << std::endl;
        std::cin.get();
    }
}

// nAlign must be a multiple of 2 and no more than 64
template <typename T, size_t nAlign>
inline T* myAlignedNew(size_t n)
{
    // allocate memory
    char *p = new char[n * sizeof(T) + nAlign];

    // align pointer
    auto addr = reinterpret_cast<std::uintptr_t>(p);
    size_t moveFwdBy = (nAlign - (addr % nAlign));
    p += moveFwdBy;
    p[-1] = (char) moveFwdBy;

    return (T*)p;
}

inline void myAlignedDelete(void *p)
{
    if (p) {
        char* p8 = (char*)p;
        p8 -= (size_t)p8[-1];
        delete [] p8;
    }
}

template <typename T, unsigned nAlign>
class AlignedVector
{
    T* m_data;
    size_t m_n;

    void deallocate()
    {
        if (m_data) {
            myAlignedDelete(m_data);
            m_data = nullptr;
        }
    }

public:
    AlignedVector() : m_data(nullptr), m_n(0) {}
    AlignedVector(size_t n) : AlignedVector() { init(n); }
    ~AlignedVector() { deallocate(); }

    AlignedVector(const AlignedVector& o) : m_data(nullptr), m_n(0)
    {
        if (o.m_n) { init(o.m_n); std::memcpy(m_data, o.m_data, m_n * sizeof(T)); }
    }
    AlignedVector(AlignedVector&& o) noexcept : m_data(o.m_data), m_n(o.m_n)
    {
        o.m_data = nullptr; o.m_n = 0;
    }
    AlignedVector& operator=(const AlignedVector& o)
    {
        if (this != &o) { init(o.m_n); std::memcpy(m_data, o.m_data, m_n * sizeof(T)); }
        return *this;
    }
    AlignedVector& operator=(AlignedVector&& o) noexcept
    {
        if (this != &o) { deallocate(); m_data = o.m_data; m_n = o.m_n; o.m_data = nullptr; o.m_n = 0; }
        return *this;
    }

    void init(size_t n)
    {
        if (n != m_n) {
            deallocate();
            m_data = myAlignedNew<T, nAlign>(n);
            m_n = n;
        }
    }

    T& operator[](size_t i) { return m_data[i]; }
    const T& operator[](size_t i) const { return m_data[i]; }
    T* data() { return m_data; }
    const T* data() const { return m_data; }
};

// Given a cube with dimension [outerDim][MidDim][InnerDim]
// Extract values with midDim=midIndex and store in a matrix of dimensions [outerDim][InnerDim]
template <size_t MidDim, size_t InnerDim, typename T>
inline void cubeToMatrix(T* vec, const T* mat, size_t outerDim, size_t midIndex)
{
    mat += InnerDim * midIndex;
    for (size_t w = 0; w < outerDim; ++w)
        for (size_t i = 0; i < InnerDim; ++i)
            vec[w * InnerDim + i] = mat[w * (MidDim * InnerDim) + i];
}

// Given a cube with dimension [outerDim][MidDim][InnerDim]
// Copy a matrix of dimensions [outerDim][InnerDim] to the layer in the cube where midDim=midIndex
template <size_t MidDim, size_t InnerDim, typename T>
inline void matrixToCube(T* mat, const T* vec, size_t outerDim, size_t midIndex)
{
    mat += InnerDim * midIndex;
    for (size_t w = 0; w < outerDim; ++w)
        for (size_t i = 0; i < InnerDim; ++i)
            mat[w * (MidDim * InnerDim) + i] = vec[w * InnerDim + i];
}

