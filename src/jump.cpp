#include "jump_matrix.h"

#include <set>
#include <thread>
#include <memory>
#include <sstream>
#include <optional>
#include <filesystem>

using namespace std;

void wait()
{
    std::cout << "press a key...";
    char c;
    std::cin >> c;
}

enum GenType { undef, mt32, mt64, sfmt };

template <GenType>
struct GenTraits;

template <>
struct GenTraits<mt32>
{
    typedef MT19937Matrix matrix_t;
    static constexpr size_t power2 = 0;
};

template <>
struct GenTraits<sfmt>
{
    typedef SFMT19937Matrix matrix_t;
    static constexpr size_t power2 = 2;
};

template <typename Matrix>
void square(const Matrix& src, Matrix& dst, std::vector<typename Matrix::buffer_t>& buffers)
{
    auto start = std::chrono::system_clock::now();
    dst.square(src, buffers);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;
}

const std::string extension = ".bits";

std::string mkFileName(const std::string& path, size_t n)
{
    std::ostringstream os;
    os << path << "F" << std::setw(5) << std::setfill('0') << n << extension;
    return os.str();
}

void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Syntax:\n"
        << "   jump -g=<generator> [-j=<nthreads>] [-p=<filepath>] [-f=<savefreq>] [-s=<stopindex>] [-wait]\n"
        << "Example:\n"
        << "   jump -g=mt32 -j=8 -p=./dat/ -f=100 -s=1000\n"
        << " generator must be one of {mt32, sfmt}\n"
        << " nthreads defaults to host concurrency\n"
        << " filepath defaults to ./dat/<genname>\n"
        << " savefrequency defaults to 100 (useful to interrupt and resume)\n"
        << " targets: comma separated list of indices to save (e.g. -t=12,52,78), defaults to {9, 100, 19933, 19934, 19935, 19936, 19937}\n";
    std::exit(-1);
}

template <GenType gen>
void run(const std::string& filepath, size_t nThreads, size_t saveFrequency, const std::set<size_t>& targets)
{
    typedef typename GenTraits<gen>::matrix_t matrix_t;
    matrix_t f[2];
    std::vector<typename matrix_t::buffer_t> buffers(nThreads);

    std::set<size_t> onDisk;
    onDisk.insert(GenTraits<gen>::power2);

    for (const auto& entry : std::filesystem::directory_iterator(filepath)) {
        if (std::filesystem::is_regular_file(entry) && entry.path().has_extension() && entry.path().extension().string() == extension) {
            std::string s = entry.path().filename().string();
            if (s.length() > extension.length() + 1 && s[0] == 'F') {
                s = s.substr(1, s.length() - 1 - extension.length());
                size_t n = (size_t)atoi(s.c_str());
                if (n >= GenTraits<gen>::power2)
                    onDisk.insert(n);
            }
        }
    }

    size_t currentIdxInMem = size_t(-1);

    for (size_t t : targets) {
        auto it = onDisk.upper_bound(t);
        if (it == onDisk.begin()) continue;
        --it;
        size_t e = *it;

        if (e == t) {
            std::cout << "Target " << t << " already exists.\n";
            continue;
        }

        if (currentIdxInMem == size_t(-1) || currentIdxInMem != e) {
            if (e == GenTraits<gen>::power2) {
                std::cout << "Initializing from base matrix: F^(2^" << e << ")\n";
                f[e % 2] = matrix_t{};
            }
            else {
                std::string fn = mkFileName(filepath, e);
                std::cout << "Loading F^(2^" << e << ") from " << fn << "\n";
                std::ifstream is(fn, std::ios::binary);
                if (is) {
                    f[e % 2].fromBin(is);
                }
                else {
                    std::cerr << "Error loading " << fn << ", falling back to base matrix\n";
                    f[e % 2] = matrix_t{};
                    e = GenTraits<gen>::power2;
                }
            }
            f[e % 2].printSparsity();
            currentIdxInMem = e;
        }

        std::string toBeDeleted = "";

        for (size_t i = e + 1; i <= t; ++i) {
            std::cout << "Computing F^(2^" << i << ")\n";
            size_t in = (i - 1) % 2;
            size_t out = i % 2;
            f[out].resetZero();
            square(f[in], f[out], buffers);
            f[out].printSparsity();
            currentIdxInMem = i;

            if ((i % saveFrequency) == 0 || i == t) {
                std::string fn = mkFileName(filepath, i);
                std::cout << "  Saving: " << fn << " ... ";
                std::ofstream of(fn, std::ios::binary);
                f[out].toBin(of);
                of.close();
                std::cout << "saved\n";
                onDisk.insert(i);

                if (!toBeDeleted.empty()) {
                    if (std::filesystem::exists(toBeDeleted)) {
                        std::filesystem::remove(toBeDeleted);
                        std::cout << "  Deleted intermediate: " << toBeDeleted << "\n";
                    }
                    toBeDeleted = "";
                }

                if (targets.find(i) == targets.end()) {
                    toBeDeleted = fn;
                }
            }
        }
    }

    std::cout << "Final cleanup: removing non-target files...\n";
    if (std::filesystem::exists(filepath)) {
        for (const auto& entry : std::filesystem::directory_iterator(filepath)) {
            if (std::filesystem::is_regular_file(entry) && entry.path().has_extension() && entry.path().extension().string() == extension) {
                std::string s = entry.path().filename().string();
                if (s.length() > extension.length() + 1 && s[0] == 'F') {
                    std::string numPart = s.substr(1, s.length() - 1 - extension.length());
                    size_t n = (size_t)atoi(numPart.c_str());
                    if (targets.find(n) == targets.end()) {
                        std::filesystem::remove(entry.path());
                        std::cout << "  Removed non-target: " << entry.path().filename().string() << "\n";
                    }
                }
            }
        }
    }
}


int main(int argc, const char** argv)
{
    ArgMap args = parseArgs(argc, argv);
    waitForDebugger(args);

    // parse command line arguments
    std::set<size_t> targetsSet;
    size_t nThreads = std::thread::hardware_concurrency();
    std::string filepath;
    std::string targets = "9,100,19933,19934,19935,19936,19937";
    std::string gentype;
    size_t saveFrequency = 100;

    try {
        // number of threads to use for matrix multiplication, defaults to host concurrency
        consumeArg(args, "j", false, nThreads);
        std::cout << "nThreads = " << nThreads << "\n";

        // generator type (mt32 or sfmt) is required
        consumeArg(args, "g", true, gentype);
        if (gentype != "mt32" && gentype != "sfmt") {
            std::cerr << "Error: invalid generator type: " << gentype << "\n";
            usage();
            return -1;
        }
        std::cout << "generator = " << gentype << "\n";

        // output directory, defaults to ./dat/<genname>/
        filepath = "./dat/" + gentype + "/";
        if (consumeArg(args, "p", false, filepath)) {
            if (filepath.back() != '/') filepath.push_back('/');
        }
        std::cout << "filepath = " << filepath << "\n";
        if (!std::filesystem::exists(filepath)) {
            std::cerr << "Error: output directory does not exist: " << filepath << "\n";
            return -1;
        }

        // save frequency (in terms of number of jumps), defaults to 100
        consumeArg(args, "f", false, saveFrequency);
        std::cout << "savefrequency = " << saveFrequency << "\n";

        // targets to save, defaults to {9, 100, 19933, 19934, 19935, 19936, 19937}
        consumeArg(args, "t", false, targets);
        {
            std::stringstream ss(targets);
            std::string index;
            while (std::getline(ss, index, ','))
                targetsSet.insert(atoi(index.c_str()));
        }
        if (!targets.empty()) {
            std::cout << "targets = ";
            for (auto it = targetsSet.begin(); it != targetsSet.end(); ++it)
                std::cout << *it << (std::next(it) == targetsSet.end() ? "" : ",");
            std::cout << "\n";
        }
        else {
            std::cerr << "Error: target indices cannot be empty\n";
            usage();
            return -1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        usage();
    }

    if (gentype == "mt32")
        run<mt32>(filepath, nThreads, saveFrequency, targetsSet);
    else if (gentype == "sfmt")
        run<sfmt>(filepath, nThreads, saveFrequency, targetsSet);
    else {
        usage();
        return -1;
    }

    return 0;

}
