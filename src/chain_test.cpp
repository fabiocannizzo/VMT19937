// chain_test — MT32/MT64/SFMT poly==matrix equivalence with in-memory and round-trip validation.
//
// At each squaring step (i=startStep+1..endStep) runs up to six tests:
//   inmem     : in-memory poly vs in-memory matrix, SIMD ISA  (no disk; catches squaring bugs)
//   scalar    : in-memory poly vs in-memory matrix, ISA::Scalar (isolates Scalar step() bugs)
//   cross_mat : in-memory poly vs disk-loaded matrix          (isolates matrix serialisation)
//   cross_poly: disk-loaded poly vs in-memory matrix          (isolates poly serialisation)
//   ident      byte-for-byte comparison against pre-existing checkpoint (NA if absent)
//   orig      : for step < 18 only: poly-jumped output vs reference generator advanced 2^n steps
//
// All output goes to stdout. Redirect / tee as needed.
//
// Thread model:
//   Main thread   — polynomial squaring + matrix squaring (N worker threads internally).
//                   Allocates a fresh Matrix per step, squares into it, then immediately
//                   pushes a WorkItem (shared_ptr<Matrix> + poly copy) to the queue and
//                   starts the next step. Never waits for I/O or tests.
//   Worker thread — persistent; pops WorkItems, saves checkpoints, runs all tests,
//                   manages rolling eviction, and logs results.
//
// Usage: chain_test [outdir] -g=mt32|mt64|sfmt [-end=N] [-start=N] [-j=N] [-keep=N] [-serial] [-h|--help]
//   outdir      directory for checkpoint files (default: dat/ subdir for chosen generator)
//   -g=type     generator type (compulsory)
//   -end=N      stop after step N (default: 19937)
//   -start=N    force start from step N (clamped to generator minimum);
//               deletes chain checkpoint files for all steps > N
//   -j=N        number of threads for matrix multiplication (default: hardware_concurrency)
//   -keep=N     number of recent chain checkpoint files to retain (default: 5)
//   -serial     enable cross_mat and cross_poly serialisation round-trip tests (default: off)
//
// Chain starts from the mathematical root: companion matrix F_0 (MT) or F_2 (SFMT) via
// default constructor, and trivial polynomial J_0=x (MT) or J_2=x (SFMT).
// All canonical F/J files in dat/matrix/ and dat/poly/ are verified via ident as the chain reaches each step.
// On success all remaining chain checkpoint files are deleted.
// On failure the checkpoint files for the failing step are preserved.

#include "RandGen.h"
#include "polynomial_jump.h"
#include "jump_matrix.h"
#include "../SFMT-src-1.5.1/SFMT.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

extern "C" {
    void               init_genrand(unsigned long s);
    unsigned long      genrand_int32();
    void               init_genrand64(unsigned long long seed);
    unsigned long long genrand64_int64();
}

// 2 full sweeps of the largest generator state (624 output words for MT32/SFMT) + 1
static const size_t g_nRand = 1249;

using namespace xvmt;
using namespace xvmt::details;
namespace fs = std::filesystem;
using namespace std::chrono;

// ── Helpers ───────────────────────────────────────────────────────────────────

static constexpr int s_colW = 5;

template <typename T>
static std::string serialize(const auto& obj, int step) {
    std::ostringstream ss;
    obj.toBin(ss, T::Gen::s_bitsGenType, (uint32_t)step);
    return ss.str();
}

static std::string readFileBytes(const std::string& path) {
    std::ifstream ifs(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(ifs), {});
}

template <typename T, typename PolyF, typename MatF>
static std::string runIdentCheck(int step, const typename T::Poly& poly, const typename T::Matrix& mat,
                                 PolyF polyFile, MatF matFile)
{
    auto cpf = T::canonPolyFile(step);
    auto cmf = T::canonMatFile(step);

    std::string poly_bytes, mat_bytes;
    bool poly_computed = false, mat_computed = false;
    auto getPolyBytes = [&]() -> const std::string& {
        if (!poly_computed) { poly_bytes = serialize<T>(poly, step); poly_computed = true; }
        return poly_bytes;
    };
    auto getMatBytes = [&]() -> const std::string& {
        if (!mat_computed) { mat_bytes = serialize<T>(mat, step); mat_computed = true; }
        return mat_bytes;
    };

    std::vector<std::string> mismatches;
    bool any_checked = false;

    if (fs::exists(polyFile(step))) {
        any_checked = true;
        if (getPolyBytes() != readFileBytes(polyFile(step))) mismatches.push_back("chain-poly");
    }
    if (fs::exists(matFile(step))) {
        any_checked = true;
        if (getMatBytes() != readFileBytes(matFile(step))) mismatches.push_back("chain-mat");
    }
    if (fs::exists(cpf)) {
        any_checked = true;
        if (getPolyBytes() != readFileBytes(cpf)) mismatches.push_back("canon-poly");
    }
    if (fs::exists(cmf)) {
        any_checked = true;
        if (getMatBytes() != readFileBytes(cmf)) mismatches.push_back("canon-mat");
    }

    if (!any_checked) return "NA";
    if (mismatches.empty()) return "OK";

    std::string r = "ERR(";
    for (size_t k = 0; k < mismatches.size(); ++k) {
        if (k > 0) r += ",";
        r += mismatches[k];
    }
    r += ")";
    return r;
}

static void printResult(const std::string& label, const std::string& res) {
    std::cout << " " << label << "=" << std::left << std::setw(s_colW) << res;
}

// ── Generator traits ──────────────────────────────────────────────────────────

struct MT32Traits {
    using Gen        = XMT19937<SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = MT19937Matrix<32>;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = XMT19937<ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint32_t;

    static constexpr int  s_reduceDegree  = 19937;
    static constexpr int  s_startStep     = 0;
    static constexpr const char* s_gentype       = "mt32";
    static constexpr const char* s_chainPrefix   = "mt32_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/mt32/characteristic.mt32.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/mt32/";

    static constexpr uint32_t s_seed = 1234UL;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/mt32/J%05d.mt32.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/mt32/F%05d.mt32.bits", step); return buf;
    }

    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static std::string runTestOrig(const PolyScalar& ps, int step) {
        init_genrand(s_seed);
        uint64_t count = uint64_t(1) << step;
        for (uint64_t k = 0; k < count; ++k) genrand_int32();
        Word ref[g_nRand];
        for (size_t k = 0; k < g_nRand; ++k) ref[k] = (uint32_t)genrand_int32();
        GenScalar gen;
        initGenPolyS(gen, ps);
        for (size_t k = 0; k < g_nRand; ++k) {
            auto a = nextWordS(gen);
            if (a != ref[k]) {
                std::ostringstream ss;
                ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k]
                   << " poly=0x" << a << std::dec << ")";
                return ss.str();
            }
        }
        return "";
    }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seed, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seed, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seed, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seed, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint32(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint32(); }
};

struct MT64Traits {
    using Gen        = XMT19937_64<SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = MT19937Matrix<64>;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = XMT19937_64<ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint64_t;

    static constexpr int  s_reduceDegree  = 19937;
    static constexpr int  s_startStep     = 0;
    static constexpr const char* s_gentype       = "mt64";
    static constexpr const char* s_chainPrefix   = "mt64_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/mt64/characteristic.mt64.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/mt64/";

    static constexpr uint64_t s_seed = 0x123456789ABCULL;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/mt64/J%05d.mt64.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/mt64/F%05d.mt64.bits", step); return buf;
    }

    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static std::string runTestOrig(const PolyScalar& ps, int step) {
        init_genrand64(s_seed);
        uint64_t count = uint64_t(1) << step;
        for (uint64_t k = 0; k < count; ++k) genrand64_int64();
        Word ref[g_nRand];
        for (size_t k = 0; k < g_nRand; ++k) ref[k] = genrand64_int64();
        GenScalar gen;
        initGenPolyS(gen, ps);
        for (size_t k = 0; k < g_nRand; ++k) {
            auto a = nextWordS(gen);
            if (a != ref[k]) {
                std::ostringstream ss;
                ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k]
                   << " poly=0x" << a << std::dec << ")";
                return ss.str();
            }
        }
        return "";
    }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seed, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seed, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seed, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seed, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint64(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint64(); }
};

struct SFMTTraits {
    using Gen        = VSFMT19937<128, false, SIMD_ISA>;
    using Poly       = Gen::poly_t;
    using PolyBig    = Polynomial<65536, SIMD_ISA>;
    using Matrix     = SFMT19937Matrix;
    using MatBuf     = Matrix::buffer_t;
    using GenScalar  = VSFMT19937<128, false, ISA::Scalar>;
    using PolyScalar = GenScalar::poly_t;
    using Word       = uint32_t;

    static constexpr int  s_reduceDegree  = 19968;
    static constexpr int  s_startStep     = 2;
    static constexpr const char* s_gentype       = "sfmt";
    static constexpr const char* s_chainPrefix   = "sfmt_chain_";
    static constexpr const char* s_charPolyFile  = "./dat/poly/sfmt/characteristic.sfmt.hex";
    static constexpr const char* s_defaultOutdir = "./dat/poly/sfmt/";

    static constexpr uint32_t s_seeds[]  = { 0x123, 0x234, 0x345, 0x456 };
    static constexpr uint32_t s_seedLen  = 4;

    static std::string canonPolyFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/poly/sfmt/J%05d.sfmt.bits", step); return buf;
    }
    static std::string canonMatFile(int step) {
        char buf[64]; snprintf(buf, sizeof(buf), "./dat/matrix/sfmt/F%05d.sfmt.bits", step); return buf;
    }

    // J_2 = x^1: one SFMT recurrence step = 4 uint32 outputs = 2^2 output words.
    static void initRootPoly(Poly& p) { p.resetZero(); p.setBit(1, true); }

    static std::string runTestOrig(const PolyScalar& ps, int step) {
        sfmt_t sfmt;
        sfmt_init_by_array(&sfmt, const_cast<uint32_t*>(s_seeds), (int)s_seedLen);
        uint64_t count = uint64_t(1) << step;
        for (uint64_t k = 0; k < count; ++k) sfmt_genrand_uint32(&sfmt);
        Word ref[g_nRand];
        for (size_t k = 0; k < g_nRand; ++k) ref[k] = sfmt_genrand_uint32(&sfmt);
        GenScalar gen;
        initGenPolyS(gen, ps);
        for (size_t k = 0; k < g_nRand; ++k) {
            auto a = nextWordS(gen);
            if (a != ref[k]) {
                std::ostringstream ss;
                ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k]
                   << " poly=0x" << a << std::dec << ")";
                return ss.str();
            }
        }
        return "";
    }

    static void initGenPoly(Gen& g, const Poly& p)              { g.reinit(s_seeds, s_seedLen, &p, nullptr); }
    static void initGenMat(Gen& g, const Matrix& m)             { g.reinit(s_seeds, s_seedLen, 1, &m, nullptr); }
    static void initGenPolyS(GenScalar& g, const PolyScalar& p) { g.reinit(s_seeds, s_seedLen, &p, nullptr); }
    static void initGenMatS(GenScalar& g, const Matrix& m)      { g.reinit(s_seeds, s_seedLen, 1, &m, nullptr); }
    static Word nextWord(Gen& g)        { return g.genrand_uint32(); }
    static Word nextWordS(GenScalar& g) { return g.genrand_uint32(); }
};

// ── WorkQueue ─────────────────────────────────────────────────────────────────

template <typename T>
struct WorkItem {
    int                                 step;
    double                              sqTime;
    typename T::Poly                    polyCopy;
    std::shared_ptr<typename T::Matrix> mat;
};

template <typename T>
class WorkQueue {
    std::queue<std::optional<WorkItem<T>>> q_;
    std::mutex                              mtx_;
    std::condition_variable                 cv_;
public:
    void push(std::optional<WorkItem<T>> item) {
        { std::lock_guard lk(mtx_); q_.push(std::move(item)); }
        cv_.notify_one();
    }
    std::optional<WorkItem<T>> pop() {
        std::unique_lock lk(mtx_);
        cv_.wait(lk, [&]{ return !q_.empty(); });
        auto v = std::move(q_.front()); q_.pop();
        return v;
    }
};

// ── Template helpers ──────────────────────────────────────────────────────────

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

template <typename T>
static std::string runTest(const typename T::Poly& poly, const typename T::Matrix& mat)
{
    typename T::Gen genPoly, genMatrix;
    T::initGenPoly(genPoly, poly);
    T::initGenMat(genMatrix, mat);
    for (size_t k = 0; k < g_nRand; ++k) {
        auto a = T::nextWord(genPoly);
        auto b = T::nextWord(genMatrix);
        if (a != b) {
            std::ostringstream ss;
            ss << "MISMATCH(word=" << k << " poly=0x" << std::hex << a
               << " mat=0x" << b << std::dec << ")";
            return ss.str();
        }
    }
    return "";
}

template <typename T>
static typename T::PolyScalar toScalarPoly(const typename T::Poly& p)
{
    alignas(64) uint32_t tmp[T::Poly::s_n32];
    p.m_data.store(tmp);
    typename T::PolyScalar ps;
    ps.m_data = typename T::PolyScalar::Reg(tmp);
    return ps;
}

template <typename T>
static std::string runTestScalar(const typename T::Poly& poly, const typename T::Matrix& mat)
{
    auto ps = toScalarPoly<T>(poly);
    typename T::GenScalar genPoly, genMatrix;
    T::initGenPolyS(genPoly, ps);
    T::initGenMatS(genMatrix, mat);
    for (size_t k = 0; k < g_nRand; ++k) {
        auto a = T::nextWordS(genPoly);
        auto b = T::nextWordS(genMatrix);
        if (a != b) {
            std::ostringstream ss;
            ss << "MISMATCH(word=" << k << " poly=0x" << std::hex << a
               << " mat=0x" << b << std::dec << ")";
            return ss.str();
        }
    }
    return "";
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

// ── runChainTest ──────────────────────────────────────────────────────────────

template <typename T>
static int runChainTest(std::string outdir, int endStep, int nThreads, int forceStart, int keepFiles, bool doSerial)
{
    if (!outdir.empty() && outdir.back() != '/') outdir += '/';

    auto polyFile = [&](int step) {
        std::ostringstream ss;
        ss << outdir << T::s_chainPrefix << "poly_"
           << std::setw(5) << std::setfill('0') << step << ".bits";
        return ss.str();
    };
    auto matFile = [&](int step) {
        std::ostringstream ss;
        ss << outdir << T::s_chainPrefix << "mat_"
           << std::setw(5) << std::setfill('0') << step << ".bits";
        return ss.str();
    };

    // Load characteristic polynomial
    typename T::Poly P;
    {
        std::ifstream ifs(T::s_charPolyFile);
        if (!ifs) { std::cerr << "Cannot open " << T::s_charPolyFile << "\n"; return 1; }
        std::string line;
        while (getline(ifs, line) && (line.empty() || line[0] == '#'));
        P.fromString(line);
    }

    std::vector<typename T::MatBuf> buffers(nThreads);

    // Scan for existing checkpoint pairs and build savedSteps
    std::deque<int> savedSteps;
    int resumeStep = -1;
    {
        std::error_code ec;
        std::vector<int> found;
        const std::string matPrefix = std::string(T::s_chainPrefix) + "mat_";
        const std::string suffix = ".bits";
        if (fs::exists(outdir)) {
            for (const auto& entry : fs::directory_iterator(outdir, ec)) {
                if (!fs::is_regular_file(entry)) continue;
                const std::string fname = entry.path().filename().string();
                if (fname.rfind(matPrefix, 0) != 0) continue;
                if (fname.size() < matPrefix.size() + suffix.size()) continue;
                if (fname.substr(fname.size() - suffix.size()) != suffix) continue;
                try {
                    int step = std::stoi(fname.substr(matPrefix.size()));
                    if (fs::exists(polyFile(step)))
                        found.push_back(step);
                } catch (...) {}
            }
        }
        std::sort(found.begin(), found.end());
        for (int s : found) {
            if (s > resumeStep) resumeStep = s;
            savedSteps.push_back(s);
        }
        while (savedSteps.size() > (size_t)keepFiles) {
            int victim = savedSteps.front(); savedSteps.pop_front();
            fs::remove(matFile(victim), ec);
            fs::remove(polyFile(victim), ec);
        }
    }

    // Determine start: forced, resumed, or fresh
    typename T::Poly poly;
    int startI = T::s_startStep;
    std::string startMode;
    auto curMat = std::make_shared<typename T::Matrix>();

    if (forceStart >= 0) {
        if (forceStart < T::s_startStep)
            forceStart = T::s_startStep;
        startI = forceStart;
        if (forceStart == T::s_startStep) {
            T::initRootPoly(poly);
            // *curMat already default-constructed
        } else {
            bool loaded = false;
            auto cpf = T::canonPolyFile(forceStart);
            auto cmf = T::canonMatFile(forceStart);
            if (fs::exists(cpf) && fs::exists(cmf)) {
                try {
                    poly = loadPoly<T>(cpf);
                    loadMat<T>(*curMat, cmf);
                    loaded = true;
                } catch (...) {}
            }
            if (!loaded) {
                try {
                    poly = loadPoly<T>(polyFile(forceStart));
                    loadMat<T>(*curMat, matFile(forceStart));
                } catch (const std::exception& e) {
                    std::cerr << "Cannot load forced start step " << forceStart << ": " << e.what() << "\n";
                    return 1;
                }
            }
        }
        startMode = "FORCED@" + std::to_string(forceStart);
    } else {
        bool resumed = (resumeStep >= T::s_startStep);
        if (resumed) {
            try {
                poly = loadPoly<T>(polyFile(resumeStep));
                loadMat<T>(*curMat, matFile(resumeStep));
                startI = resumeStep;
            } catch (...) {
                resumed = false;
            }
        }
        if (!resumed) {
            T::initRootPoly(poly);
            // *curMat already default-constructed
            startI = T::s_startStep;
        }
        startMode = resumed ? "RESUMED" : "FRESH";
    }

    // When forced: remove chain checkpoint files for all steps after the start point
    if (forceStart >= 0) {
        std::error_code ec2;
        savedSteps.erase(
            std::remove_if(savedSteps.begin(), savedSteps.end(),
                [&](int s) {
                    if (s > forceStart) {
                        fs::remove(matFile(s), ec2);
                        fs::remove(polyFile(s), ec2);
                        return true;
                    }
                    return false;
                }),
            savedSteps.end());
    }

    std::cout << "chain_test: " << startMode
              << "  g=" << T::s_gentype
              << "  startI=" << startI
              << "  endI=" << endStep
              << "  threads=" << nThreads
              << "  g_nRand=" << g_nRand
              << "  keep=" << keepFiles << "\n"
              << "#\n"
              << "# Output columns per step:\n"
              << "#   i          squaring index\n"
              << "#   inmem      SIMD poly and matrix both kept in RAM\n"
              << "#   scalar     same as inmem but with ISA::Scalar generator\n"
              << "#   cross_mat  in-memory poly vs disk-loaded matrix       (NA unless -serial)\n"
              << "#   cross_poly disk-loaded poly vs in-memory matrix       (NA unless -serial)\n"
              << "#   ident      byte-for-byte match against pre-existing checkpoint (NA if absent)\n"
              << "#   orig       step < 18 only: poly-jumped output vs reference generator run for 2^n steps\n"
              << "#   t          matrix + poly squaring time (main thread)\n"
              << "#   eta        estimated time remaining (based on this session)\n"
              << "#\n" << std::flush;

    // Validate starting point synchronously
    {
        auto t0 = steady_clock::now();
        auto r_inmem  = runTest<T>(poly, *curMat);
        auto r_scalar = runTestScalar<T>(poly, *curMat);
        auto r_ident  = runIdentCheck<T>(startI, poly, *curMat, polyFile, matFile);
        std::string r_orig = "NA";
        if (startI < 18) {
            auto ps = toScalarPoly<T>(poly);
            r_orig = T::runTestOrig(ps, startI);
        }
        if (r_orig.empty()) r_orig = "OK";

        double dt = duration<double>(steady_clock::now() - t0).count();
        std::cout << "i=" << std::right << std::setw(5) << std::setfill('0') << startI << std::setfill(' ');
        printResult("inmem",      r_inmem.empty()  ? "OK" : r_inmem);
        printResult("scalar",     r_scalar.empty() ? "OK" : r_scalar);
        printResult("cross_mat",  "NA");
        printResult("cross_poly", "NA");
        printResult("ident",      r_ident);
        printResult("orig",       r_orig);
        std::cout << " t=" << std::fixed << std::setprecision(3) << dt << "s\n" << std::flush;

        if (!r_inmem.empty() || !r_scalar.empty() || (r_orig != "OK" && r_orig != "NA")) {
            std::cout << "STOPPING: starting point fails\n";
            return 2;
        }
    }

    std::atomic<bool> workerFailed{false};
    auto sessionStart = steady_clock::now();
    std::atomic<int>  sessionIters{0};

    WorkQueue<T> workQueue;

    // Worker thread: persistent, never returns until poison pill
    std::thread worker([&]{
        std::error_code ec;
        for (;;) {
            auto opt = workQueue.pop();
            if (!opt) break;
            WorkItem<T> w = std::move(*opt);

            // Ident: byte-for-byte check against any pre-existing reference files.
            std::string r_ident = runIdentCheck<T>(w.step, w.polyCopy, *w.mat, polyFile, matFile);

            if (!fs::exists(outdir)) fs::create_directories(outdir, ec);
            saveMat<T>(*w.mat,      matFile(w.step), w.step);
            savePoly<T>(w.polyCopy, polyFile(w.step), w.step);

            auto r_inmem      = runTest<T>(w.polyCopy, *w.mat);
            auto r_scalar     = runTestScalar<T>(w.polyCopy, *w.mat);

            std::string r_cross_mat = "NA", r_cross_poly = "NA";
            if (doSerial) {
                typename T::Poly   poly_rt = loadPoly<T>(polyFile(w.step));
                typename T::Matrix mat_rt;  loadMat<T>(mat_rt, matFile(w.step));
                r_cross_mat  = runTest<T>(w.polyCopy, mat_rt);
                r_cross_poly = runTest<T>(poly_rt,   *w.mat);
            }

            std::string r_orig;
            if (w.step < 18) {
                auto ps = toScalarPoly<T>(w.polyCopy);
                r_orig = T::runTestOrig(ps, w.step);
            }

            bool ok = r_inmem.empty() && r_scalar.empty() &&
                      (r_cross_mat.empty()  || r_cross_mat  == "NA") &&
                      (r_cross_poly.empty() || r_cross_poly == "NA") &&
                      (r_ident == "OK" || r_ident == "NA") &&
                      r_orig.empty();

            if (ok) {
                savedSteps.push_back(w.step);
                if (savedSteps.size() > (size_t)keepFiles) {
                    int victim = savedSteps.front(); savedSteps.pop_front();
                    fs::remove(matFile(victim), ec);
                    fs::remove(polyFile(victim), ec);
                }
            }
            int iters = ++sessionIters;
            double elapsed = duration<double>(steady_clock::now() - sessionStart).count();

            std::cout << "i=" << std::right << std::setw(5) << std::setfill('0') << w.step << std::setfill(' ');
            printResult("inmem",      r_inmem.empty()      ? "OK" : r_inmem);
            printResult("scalar",     r_scalar.empty()     ? "OK" : r_scalar);
            printResult("cross_mat",  r_cross_mat.empty()  ? "OK" : r_cross_mat);
            printResult("cross_poly", r_cross_poly.empty() ? "OK" : r_cross_poly);
            printResult("ident",      r_ident);

            std::string r_orig_disp = "NA";
            if (w.step < 18) {
                r_orig_disp = r_orig.empty() ? "OK" : r_orig;
            }
            printResult("orig", r_orig_disp);

            std::cout << " t=" << std::fixed << std::setprecision(3) << w.sqTime << "s";
            if (iters > 0) {
                double etaSecs = elapsed / iters * (endStep - w.step);
                int eta_h = (int)(etaSecs / 3600);
                int eta_m = (int)(etaSecs / 60) % 60;
                int eta_s = (int)etaSecs % 60;
                std::cout << "  eta=" << eta_h << "h"
                          << std::setw(2) << std::setfill('0') << eta_m << "m"
                          << std::setw(2) << std::setfill('0') << eta_s << "s";
            }
            std::cout << "\n" << std::flush;

            if (!ok) {
                workerFailed.store(true);
                break;
            }
        }
    });

    // Main squaring loop
    for (int i = startI + 1; i <= endStep && !workerFailed.load(std::memory_order_relaxed); ++i) {
        auto matNew = std::make_shared<typename T::Matrix>();

        auto tSq0 = steady_clock::now();
        poly = poly_sq_mod<T>(poly, P);
        matNew->square(*curMat, buffers);
        double sqTime = duration<double>(steady_clock::now() - tSq0).count();

        workQueue.push(WorkItem<T>{i, sqTime, poly, matNew});
        curMat = std::move(matNew);
    }

    workQueue.push(std::nullopt);
    worker.join();

    if (!workerFailed.load()) {
        std::cout << "ALL OK: no mismatch for i=" << startI << ".." << endStep << "\n" << std::flush;
        std::error_code ec;
        if (fs::exists(outdir)) {
            for (const auto& entry : fs::directory_iterator(outdir, ec)) {
                if (!fs::is_regular_file(entry)) continue;
                if (entry.path().filename().string().rfind(T::s_chainPrefix, 0) == 0)
                    fs::remove(entry.path(), ec);
            }
        }
    } else {
        std::cout << "STOPPED due to failure\n" << std::flush;
    }

    return workerFailed.load() ? 3 : 0;
}

// ── main ──────────────────────────────────────────────────────────────────────

static void usage(const char* prog)
{
    unsigned hw = std::thread::hardware_concurrency();
    std::cerr <<
        "Usage: " << prog << " [outdir] -g=mt32|mt64|sfmt [-end=N] [-start=N] [-j=N] [-keep=N] [-serial] [-h|--help]\n"
        "\n"
        "  outdir          directory for checkpoint files\n"
        "                    default: ./dat/poly/mt32/ or ./dat/poly/mt64/ or ./dat/poly/sfmt/\n"
        "  -g=mt32         run MT19937-32 chain (starts from step 0)\n"
        "  -g=mt64         run MT19937-64 chain (starts from step 0)\n"
        "  -g=sfmt         run SFMT19937 chain (starts from step 2)\n"
        "  -end=N          stop after squaring step N (default: 19937)\n"
        "  -start=N        force start from step N (clamped to min for generator); deletes chain files for steps > N\n"
        "  -j=N            threads for matrix multiplication (default: " << (hw ? hw : 1) << ")\n"
        "  -keep=N         recent chain checkpoint files to retain (default: 5)\n"
        "  -serial         enable cross_mat and cross_poly serialisation round-trip tests (default: off)\n"
        "  -h, --help      show this message\n"
        "\n"
        "Chain starts from root: F_0/J_0=x for MT, F_2/J_2=x for SFMT.\n"
        "All canonical F/J files in dat/matrix/ and dat/poly/ are verified via ident as the chain passes them.\n"
        "Tests per step:\n"
        "  inmem      SIMD poly vs matrix, both in RAM\n"
        "  scalar     same with ISA::Scalar generator\n"
        "  cross_mat  in-memory poly vs disk-loaded matrix  (NA unless -serial)\n"
        "  cross_poly disk-loaded poly vs in-memory matrix  (NA unless -serial)\n"
        "  ident      byte-for-byte match vs chain checkpoint and/or canonical F/J files (NA if absent)\n"
        "  orig       step < 18: poly-jumped output vs reference generator run for 2^n steps\n"
        "On clean completion all remaining chain checkpoint files are deleted.\n";
}

int main(int argc, const char** argv)
{
    std::string gentype;
    int endStep    = 19937;
    int forceStart = -1;
    int nThreads   = (int)std::max(1u, std::thread::hardware_concurrency());
    int keepFiles  = 5;
    bool doSerial  = false;
    std::string outdir;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "-h" || arg == "--help")    { usage(argv[0]); return 0; }
        else if (arg.rfind("-g=",       0) == 0)    gentype    = arg.substr(3);
        else if (arg.rfind("-end=",     0) == 0)    endStep    = std::stoi(arg.substr(5));
        else if (arg.rfind("-start=",   0) == 0)    forceStart = std::stoi(arg.substr(7));
        else if (arg.rfind("-j=",       0) == 0)    nThreads   = std::stoi(arg.substr(3));
        else if (arg.rfind("-threads=", 0) == 0)    nThreads   = std::stoi(arg.substr(9));
        else if (arg.rfind("-keep=",    0) == 0)    keepFiles  = std::stoi(arg.substr(6));
        else if (arg == "-serial")                  doSerial   = true;
        else if (arg[0] == '-') { std::cerr << "Unknown option: " << arg << "\n\n"; usage(argv[0]); return 1; }
        else outdir = arg;
    }

    if (gentype.empty()) {
        std::cerr << "Error: Generator type (-g) is required.\n\n";
        usage(argv[0]);
        return 1;
    }

    if (gentype == "mt32") {
        if (outdir.empty()) outdir = MT32Traits::s_defaultOutdir;
        return runChainTest<MT32Traits>(outdir, endStep, nThreads, forceStart, keepFiles, doSerial);
    } else if (gentype == "mt64") {
        if (outdir.empty()) outdir = MT64Traits::s_defaultOutdir;
        return runChainTest<MT64Traits>(outdir, endStep, nThreads, forceStart, keepFiles, doSerial);
    } else if (gentype == "sfmt") {
        if (outdir.empty()) outdir = SFMTTraits::s_defaultOutdir;
        return runChainTest<SFMTTraits>(outdir, endStep, nThreads, forceStart, keepFiles, doSerial);
    } else {
        std::cerr << "Error: Unknown generator type: " << gentype << "\n\n";
        usage(argv[0]);
        return 1;
    }
}
