#ifndef XVMT_JUMP_UTILS_H
#define XVMT_JUMP_UTILS_H

#include "jump_traits.h"
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

namespace xvmt {
namespace details {

template <typename T>
static typename T::Poly poly_sq_mod(const typename T::Poly& in, const typename T::Poly& P)
{
    typename T::PolyBig sq;
    PolyOps::square(in, sq);
    PolyOps::reduce<T::s_reduceDegree>(sq, P);
    alignas(64) uint32_t tmp[65536 / 32];
    sq.m_data.store(tmp);
    typename T::Poly res;
    res.m_data = typename T::Poly::Reg(tmp);
    return res;
}

static inline std::string readFileBytes(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(ifs), {});
}

template <typename T>
static void savePoly(const typename T::Poly& poly, const std::string& path, int step)
{
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot write " + path);
    poly.toBin(ofs, T::Gen::s_bitsGenType, (uint32_t)step);
}

template <typename T>
static typename T::Poly loadPoly(const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot read " + path);
    typename T::Poly p; p.fromBinStream(ifs);
    return p;
}

template <typename T>
static void saveMat(const typename T::Matrix& mat, const std::string& path, int step)
{
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot write " + path);
    mat.toBin(ofs, T::Gen::s_bitsGenType, (uint32_t)step);
}

template <typename T>
static void loadMat(typename T::Matrix& mat, const std::string& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot read " + path);
    mat.fromBinStream(ifs);
}

template <typename T>
static std::string poly_to_string(const typename T::Poly& p) {
    return p.toString();
}

} // namespace details
} // namespace xvmt

#endif // XVMT_JUMP_UTILS_H
