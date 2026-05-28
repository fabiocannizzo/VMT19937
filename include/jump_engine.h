#ifndef XVMT_JUMP_ENGINE_H
#define XVMT_JUMP_ENGINE_H

#include "jump_traits.h"
#include "jump_utils.h"
#include <functional>
#include <memory>
#include <chrono>
#include <filesystem>
#include <algorithm>
#include <vector>

namespace xvmt {
namespace details {

namespace fs = std::filesystem;
using namespace std::chrono;

template <typename T>
class JumpEngine {
public:
    using Poly = typename T::Poly;
    using Matrix = typename T::Matrix;
    using MatBuf = typename T::MatBuf;
    using MatrixPtr = std::shared_ptr<Matrix>;

    using OnStepCallback = std::function<void(int step, const Poly&, const MatrixPtr&, double sqTime)>;

    struct Config {
        std::string outdir;
        std::string charPolyFile;
        int         nThreads;
        int         forceStart = -1;
        std::string chainPrefix; // Optional: for chain_test style checkpoints
    };

    static int run(const Config& cfg, int endStep, OnStepCallback onStep) {
        // Load characteristic polynomial
        Poly P;
        {
            std::ifstream ifs(cfg.charPolyFile.empty() ? T::s_charPolyFile : cfg.charPolyFile);
            if (!ifs) { 
                std::cerr << "Cannot open characteristic polynomial file: " 
                          << (cfg.charPolyFile.empty() ? T::s_charPolyFile : cfg.charPolyFile) << "\n"; 
                return 1; 
            }
            std::string line;
            while (getline(ifs, line) && (line.empty() || line[0] == '#'));
            P.fromString(line);
        }

        std::vector<MatBuf> buffers(cfg.nThreads);

        Poly poly;
        MatrixPtr curMat = std::make_shared<Matrix>();
        int startI = T::s_startStep;

        // Resume logic
        if (cfg.forceStart >= 0) {
            startI = std::max(T::s_startStep, cfg.forceStart);
            if (startI == T::s_startStep) {
                T::initRootPoly(poly);
            } else {
                if (!loadStep(cfg, startI, poly, *curMat)) {
                    std::cerr << "Cannot load forced start step " << startI << "\n";
                    return 1;
                }
            }
        } else {
            int resumeStep = findHighestStep(cfg);
            if (resumeStep >= T::s_startStep) {
                if (loadStep(cfg, resumeStep, poly, *curMat)) {
                    startI = resumeStep;
                } else {
                    T::initRootPoly(poly);
                    startI = T::s_startStep;
                }
            } else {
                T::initRootPoly(poly);
                startI = T::s_startStep;
            }
        }

        // Notify caller of the starting state
        onStep(startI, poly, curMat, 0.0);

        // Main loop
        for (int i = startI + 1; i <= endStep; ++i) {
            auto matNew = std::make_shared<Matrix>();

            auto tSq0 = steady_clock::now();
            poly = poly_sq_mod<T>(poly, P);
            matNew->square(*curMat, buffers);
            double sqTime = duration<double>(steady_clock::now() - tSq0).count();

            onStep(i, poly, matNew, sqTime);
            curMat = std::move(matNew);
        }

        return 0;
    }

private:
    static std::string getPolyPath(const Config& cfg, int step) {
        if (!cfg.chainPrefix.empty()) {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s%spoly_%05d.bits", cfg.outdir.c_str(), cfg.chainPrefix.c_str(), step);
            return buf;
        }
        char buf[128];
        snprintf(buf, sizeof(buf), "%sJ%05d.%s.bits", cfg.outdir.c_str(), step, T::s_gentype);
        return buf;
    }

    static std::string getMatPath(const Config& cfg, int step) {
        if (!cfg.chainPrefix.empty()) {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s%smat_%05d.bits", cfg.outdir.c_str(), cfg.chainPrefix.c_str(), step);
            return buf;
        }
        char buf[128];
        snprintf(buf, sizeof(buf), "%sF%05d.%s.bits", cfg.outdir.c_str(), step, T::s_gentype);
        return buf;
    }

    static bool loadStep(const Config& cfg, int step, Poly& p, Matrix& m) {
        try {
            std::string pp = getPolyPath(cfg, step);
            std::string mp = getMatPath(cfg, step);
            
            // Try config paths first
            if (fs::exists(pp) && fs::exists(mp)) {
                p = loadPoly<T>(pp);
                loadMat<T>(m, mp);
                return true;
            }
            
            // Fallback to canonical paths if not found in outdir
            std::string cpp = T::canonPolyFile(step);
            std::string cmp = T::canonMatFile(step);
            if (fs::exists(cpp) && fs::exists(cmp)) {
                p = loadPoly<T>(cpp);
                loadMat<T>(m, cmp);
                return true;
            }
        } catch (...) {}
        return false;
    }

    static int findHighestStep(const Config& cfg) {
        int highest = -1;
        if (!fs::exists(cfg.outdir)) return -1;
        
        for (const auto& entry : fs::directory_iterator(cfg.outdir)) {
            if (!entry.is_regular_file()) continue;
            std::string fname = entry.path().filename().string();
            
            int step = -1;
            if (!cfg.chainPrefix.empty()) {
                std::string matPrefix = cfg.chainPrefix + "mat_";
                if (fname.rfind(matPrefix, 0) == 0) {
                    try { step = std::stoi(fname.substr(matPrefix.size())); } catch (...) {}
                }
            } else {
                if (fname.size() >= 7 && fname[0] == 'F' && fname.find(T::s_gentype) != std::string::npos) {
                    try { step = std::stoi(fname.substr(1, 5)); } catch (...) {}
                }
            }
            
            if (step >= 0) {
                // Verify both poly and mat exist
                if (fs::exists(getPolyPath(cfg, step)) && fs::exists(getMatPath(cfg, step))) {
                    if (step > highest) highest = step;
                }
            }
        }
        return highest;
    }
};

} // namespace details
} // namespace xvmt

#endif // XVMT_JUMP_ENGINE_H
