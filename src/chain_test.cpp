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
//   Main thread   — polynomial squaring + matrix squaring (via JumpEngine).
//                   Pushes a WorkItem (shared_ptr<Matrix> + poly copy) to the queue and
//                   starts the next step. Never waits for I/O or tests.
//   Worker thread — persistent; pops WorkItems, saves checkpoints, runs all tests,
//                   manages rolling eviction, and logs results.

#include "jump_engine.h"
#include "cli_args.h"

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

// 2 full sweeps of the largest generator state (624 output words for MT32/SFMT) + 1
static const size_t g_nRand = 1249;

using namespace xvmt;
using namespace xvmt::details;
namespace fs = std::filesystem;
using namespace std::chrono;

static constexpr int s_colW = 5;

// ── Helpers ───────────────────────────────────────────────────────────────────

template <typename T>
static std::string serialize(const auto& obj, int step) {
    std::ostringstream ss;
    obj.toBin(ss, T::Gen::s_bitsGenType, (uint32_t)step);
    return ss.str();
}

template <typename T>
static std::string runIdentCheck(int step, const typename T::Poly& poly, const typename T::Matrix& mat,
                                 const std::string& outdir)
{
    auto polyFile = [&](int s) {
        char buf[128]; snprintf(buf, sizeof(buf), "%s%spoly_%05d.bits", outdir.c_str(), T::s_chainPrefix, s); return std::string(buf);
    };
    auto matFile = [&](int s) {
        char buf[128]; snprintf(buf, sizeof(buf), "%s%smat_%05d.bits", outdir.c_str(), T::s_chainPrefix, s); return std::string(buf);
    };

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

static std::string runTestOrig_MT32(const MT32Traits::PolyScalar& ps, int step) {
    init_genrand(MT32Traits::s_seed);
    uint64_t count = uint64_t(1) << step;
    for (uint64_t k = 0; k < count; ++k) genrand_int32();
    uint32_t ref[g_nRand];
    for (size_t k = 0; k < g_nRand; ++k) ref[k] = (uint32_t)genrand_int32();
    MT32Traits::GenScalar gen;
    MT32Traits::initGenPolyS(gen, ps);
    for (size_t k = 0; k < g_nRand; ++k) {
        auto a = gen.genrand_uint32();
        if (a != ref[k]) {
            std::ostringstream ss;
            ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k] << " poly=0x" << a << std::dec << ")";
            return ss.str();
        }
    }
    return "";
}

static std::string runTestOrig_MT64(const MT64Traits::PolyScalar& ps, int step) {
    init_genrand64(MT64Traits::s_seed);
    uint64_t count = uint64_t(1) << step;
    for (uint64_t k = 0; k < count; ++k) genrand64_int64();
    uint64_t ref[g_nRand];
    for (size_t k = 0; k < g_nRand; ++k) ref[k] = genrand64_int64();
    MT64Traits::GenScalar gen;
    MT64Traits::initGenPolyS(gen, ps);
    for (size_t k = 0; k < g_nRand; ++k) {
        auto a = gen.genrand_uint64();
        if (a != ref[k]) {
            std::ostringstream ss;
            ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k] << " poly=0x" << a << std::dec << ")";
            return ss.str();
        }
    }
    return "";
}

static std::string runTestOrig_SFMT(const SFMTTraits::PolyScalar& ps, int step) {
    sfmt_t sfmt;
    sfmt_init_by_array(&sfmt, const_cast<uint32_t*>(SFMTTraits::s_seeds), (int)SFMTTraits::s_seedLen);
    uint64_t count = uint64_t(1) << step;
    for (uint64_t k = 0; k < count; ++k) sfmt_genrand_uint32(&sfmt);
    uint32_t ref[g_nRand];
    for (size_t k = 0; k < g_nRand; ++k) ref[k] = sfmt_genrand_uint32(&sfmt);
    SFMTTraits::GenScalar gen;
    SFMTTraits::initGenPolyS(gen, ps);
    for (size_t k = 0; k < g_nRand; ++k) {
        auto a = gen.genrand_uint32();
        if (a != ref[k]) {
            std::ostringstream ss;
            ss << "MISMATCH(word=" << k << " orig=0x" << std::hex << ref[k] << " poly=0x" << a << std::dec << ")";
            return ss.str();
        }
    }
    return "";
}

template <typename T>
static std::string runTestOrig(const typename T::PolyScalar& ps, int step) {
    if constexpr (std::is_same_v<T, MT32Traits>) return runTestOrig_MT32(ps, step);
    else if constexpr (std::is_same_v<T, MT64Traits>) return runTestOrig_MT64(ps, step);
    else if constexpr (std::is_same_v<T, SFMTTraits>) return runTestOrig_SFMT(ps, step);
    return "NA";
}

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

// ── runChainTest ──────────────────────────────────────────────────────────────

template <typename T>
static int runChainTest(std::string outdir, int endStep, int nThreads, int forceStart, int keepFiles, bool doSerial)
{
    if (!outdir.empty() && outdir.back() != '/' && outdir.back() != '\\') outdir += '/';

    auto polyFile = [&](int s) {
        char buf[128]; snprintf(buf, sizeof(buf), "%s%spoly_%05d.bits", outdir.c_str(), T::s_chainPrefix, s); return std::string(buf);
    };
    auto matFile = [&](int s) {
        char buf[128]; snprintf(buf, sizeof(buf), "%s%smat_%05d.bits", outdir.c_str(), T::s_chainPrefix, s); return std::string(buf);
    };

    std::cout << "chain_test: " 
              << "  g=" << T::s_gentype
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

    std::atomic<bool> workerFailed{false};
    auto sessionStart = steady_clock::now();
    std::atomic<int>  sessionIters{0};
    std::deque<int>   savedSteps;
    std::mutex        savedMtx;

    WorkQueue<T> workQueue;

    // Worker thread: persistent, never returns until poison pill
    std::thread worker([&]{
        std::error_code ec;
        for (;;) {
            auto opt = workQueue.pop();
            if (!opt) break;
            WorkItem<T> w = std::move(*opt);

            if (w.sqTime == 0.0) { // Starting step notification
                auto r_inmem  = runTest<T>(w.polyCopy, *w.mat);
                auto r_scalar = runTestScalar<T>(w.polyCopy, *w.mat);
                auto r_ident  = runIdentCheck<T>(w.step, w.polyCopy, *w.mat, outdir);
                std::string r_orig = "NA";
                if (w.step < 18) {
                    auto ps = toScalarPoly<T>(w.polyCopy);
                    r_orig = runTestOrig<T>(ps, w.step);
                }
                if (r_orig.empty()) r_orig = "OK";

                std::cout << "i=" << std::right << std::setw(5) << std::setfill('0') << w.step << std::setfill(' ');
                printResult("inmem",      r_inmem.empty()  ? "OK" : r_inmem);
                printResult("scalar",     r_scalar.empty() ? "OK" : r_scalar);
                printResult("cross_mat",  "NA");
                printResult("cross_poly", "NA");
                printResult("ident",      r_ident);
                printResult("orig",       r_orig);
                std::cout << " t=0.000s\n" << std::flush;

                if (!r_inmem.empty() || !r_scalar.empty() || (r_orig != "OK" && r_orig != "NA")) {
                    workerFailed.store(true); break;
                }
                continue;
            }

            // Ident: byte-for-byte check against any pre-existing reference files.
            std::string r_ident = runIdentCheck<T>(w.step, w.polyCopy, *w.mat, outdir);

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
                r_orig = runTestOrig<T>(ps, w.step);
            }

            bool ok = r_inmem.empty() && r_scalar.empty() &&
                      (r_cross_mat.empty()  || r_cross_mat  == "NA") &&
                      (r_cross_poly.empty() || r_cross_poly == "NA") &&
                      (r_ident == "OK" || r_ident == "NA") &&
                      r_orig.empty();

            if (ok) {
                std::lock_guard lk(savedMtx);
                savedSteps.push_back(w.step);
                while (savedSteps.size() > (size_t)keepFiles) {
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
                std::cout << "  eta=" << eta_h << "h" << std::setw(2) << std::setfill('0') << eta_m << "m" << std::setw(2) << std::setfill('0') << eta_s << "s";
            }
            std::cout << "\n" << std::flush;

            if (!ok) {
                workerFailed.store(true);
                break;
            }
        }
    });

    // Invoke JumpEngine
    typename JumpEngine<T>::Config cfg;
    cfg.outdir = outdir;
    cfg.nThreads = nThreads;
    cfg.forceStart = forceStart;
    cfg.chainPrefix = T::s_chainPrefix;

    auto onStep = [&](int step, const typename T::Poly& poly, const typename JumpEngine<T>::MatrixPtr& mat, double sqTime) {
        workQueue.push(WorkItem<T>{step, sqTime, poly, mat});
    };

    JumpEngine<T>::run(cfg, endStep, onStep);

    workQueue.push(std::nullopt);
    worker.join();

    if (!workerFailed.load()) {
        std::cout << "ALL OK: chain validated up to i=" << endStep << "\n" << std::flush;
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
        "  -start=N        force start from step N (clamped to min for generator)\n"
        "  -j=N            threads for matrix multiplication (default: " << (hw ? hw : 1) << ")\n"
        "  -keep=N         recent chain checkpoint files to retain (default: 5)\n"
        "  -serial         enable cross_mat and cross_poly serialisation round-trip tests (default: off)\n"
        "  -h, --help      show this message\n";
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
