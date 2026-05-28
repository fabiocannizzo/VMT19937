#include "jump_engine.h"
#include "cli_args.h"
#include <set>
#include <iostream>
#include <sstream>

using namespace xvmt;
using namespace xvmt::details;

enum Mode { mode_both, mode_matrix, mode_poly };

template <typename T>
void runGenerator(const std::string& gentype, const typename JumpEngine<T>::Config& cfg, int endStep, 
                  Mode mode, const std::set<int>& targets, int saveFrequency) 
{
    auto startTime = steady_clock::now();
    int nComputed = 0;

    auto onStep = [&](int step, const typename T::Poly& poly, const typename JumpEngine<T>::MatrixPtr& mat, double sqTime) {
        bool isTarget = targets.count(step);
        bool isFreq = (saveFrequency > 0 && (step % saveFrequency) == 0);
        bool isEnd = (step == endStep);

        if (step > T::s_startStep) {
            nComputed++;
            double elapsed = duration<double>(steady_clock::now() - startTime).count();
            std::cout << "i=" << std::setw(5) << std::setfill('0') << step 
                      << " sq=" << std::fixed << std::setprecision(2) << sqTime << "s"
                      << " elapsed=" << std::fixed << std::setprecision(1) << elapsed << "s";
            if (nComputed > 0) {
                double eta = elapsed / nComputed * (endStep - step);
                std::cout << " eta=" << (int)(eta / 3600) << "h" << std::setw(2) << std::setfill('0') << (int)(eta / 60) % 60 << "m";
            }
            if (isTarget) std::cout << " [TARGET]";
            std::cout << "\n";
        }

        if (isTarget || isFreq || isEnd) {
            std::ostringstream osp, osm;
            osp << cfg.outdir << "J" << std::setw(5) << std::setfill('0') << step << "." << T::s_gentype << ".bits";
            osm << cfg.outdir << "F" << std::setw(5) << std::setfill('0') << step << "." << T::s_gentype << ".bits";
            std::string pp = osp.str();
            std::string mp = osm.str();
            
            if (mode == mode_both || mode == mode_poly) {
                savePoly<T>(poly, pp, step);
                // std::cout << "  Saved poly: " << pp << "\n";
            }
            if (mode == mode_both || mode == mode_matrix) {
                saveMat<T>(*mat, mp, step);
                // std::cout << "  Saved mat:  " << mp << "\n";
            }
        }
    };

    JumpEngine<T>::run(cfg, endStep, onStep);
}

void usage() {
    std::cerr << "Usage: jump_generator -g=<mt32|mt64|sfmt> [-mode=<both|matrix|poly>] [-t=<targets>] [-end=<N>] [-f=<freq>] [-j=<threads>] [-outdir=<dir>] [-charpoly=<file>]\n"
              << "  -g: generator type (required)\n"
              << "  -mode: output mode (default: both)\n"
              << "  -t: comma-separated list of target steps (e.g. 19937,9)\n"
              << "  -end: stop after step N (default: highest target or 19937)\n"
              << "  -f: save frequency (default: 0, only save targets and end)\n"
              << "  -j: number of threads (default: hardware concurrency)\n"
              << "  -outdir: output directory (default: ./dat/poly/<gentype>/)\n"
              << "  -charpoly: characteristic polynomial file\n";
    std::exit(1);
}

int main(int argc, const char** argv) {
    ArgMap args = parseArgs(argc, argv);
    if (consumeArg(args, "h") || consumeArg(args, "help")) usage();

    std::string gentype;
    if (!consumeArg(args, "g", true, gentype)) usage();

    std::string modeStr = "both";
    consumeArg(args, "mode", false, modeStr);
    Mode mode = mode_both;
    if (modeStr == "matrix") mode = mode_matrix;
    else if (modeStr == "poly") mode = mode_poly;

    std::string targetsStr;
    std::set<int> targets;
    if (consumeArg(args, "t", false, targetsStr)) {
        std::stringstream ss(targetsStr);
        std::string item;
        while (std::getline(ss, item, ',')) {
            try { targets.insert(std::stoi(item)); } catch (...) {}
        }
    }

    int endStep = 19937;
    if (!targets.empty()) endStep = *targets.rbegin();
    consumeArg(args, "end", false, endStep);

    int saveFrequency = 0;
    consumeArg(args, "f", false, saveFrequency);

    int nThreads = (int)std::thread::hardware_concurrency();
    consumeArg(args, "j", false, nThreads);

    std::string outdir = "./dat/poly/" + gentype + "/";
    if (consumeArg(args, "outdir", false, outdir)) {
        if (!outdir.empty() && outdir.back() != '/' && outdir.back() != '\\') outdir += '/';
    }

    std::string charPoly;
    consumeArg(args, "charpoly", false, charPoly);

    std::cout << "jump_generator: g=" << gentype << " mode=" << modeStr << " end=" << endStep << " threads=" << nThreads << "\n";

    if (gentype == "mt32") {
        typename JumpEngine<MT32Traits>::Config cfg;
        cfg.nThreads = nThreads;
        cfg.outdir = outdir;
        cfg.charPolyFile = charPoly;
        runGenerator<MT32Traits>(gentype, cfg, endStep, mode, targets, saveFrequency);
    } else if (gentype == "mt64") {
        typename JumpEngine<MT64Traits>::Config cfg;
        cfg.nThreads = nThreads;
        cfg.outdir = outdir;
        cfg.charPolyFile = charPoly;
        runGenerator<MT64Traits>(gentype, cfg, endStep, mode, targets, saveFrequency);
    } else if (gentype == "sfmt") {
        typename JumpEngine<SFMTTraits>::Config cfg;
        cfg.nThreads = nThreads;
        cfg.outdir = outdir;
        cfg.charPolyFile = charPoly;
        runGenerator<SFMTTraits>(gentype, cfg, endStep, mode, targets, saveFrequency);
    } else { std::cerr << "Unknown generator: " << gentype << "\n"; return 1; }

    return 0;
}
