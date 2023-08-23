#include "jump_matrix.h"

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

enum GenType { undef, mt, sfmt };

template <GenType>
struct GenTraits;

template <>
struct GenTraits<mt>
{
    typedef MT19937Matrix matrix_t;
    static const size_t power2 = 0;
};

template <>
struct GenTraits<sfmt>
{
    typedef SFMT19937Matrix matrix_t;
    static const size_t power2 = 2;
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
        << "Example:\n"
        << "jump -g generator [-j nthreads] [-p filepath] [-f savefrequency] [-s stopindex]\n"
        << "  nthreads defaults to 1\n"
        << "  filepath defaults to ./\n"
        << "  savefrequency defaults to 100\n"
        << "  generator must be one of {mt, sfmt}"
        << "  stopindex: stops after saving the matrix with power >= stopIndex\n";
    std::exit(-1);
}

template <GenType gen>
void run(const std::string& filepath, size_t nThreads, size_t saveFrequency, size_t stopIndex)
{
    typename GenTraits<gen>::matrix_t f[2];

    size_t lastComputed = GenTraits<gen>::power2;

    for (const auto& entry : std::filesystem::directory_iterator(filepath)) {
        if (std::filesystem::is_regular_file(entry) && entry.path().has_extension() && entry.path().extension().string() == extension) {
            std::string s = entry.path().filename().string();
            s = s.substr(1, s.length() - 1 - extension.length());
            int n = atoi(s.c_str());
            if (n > lastComputed)
                lastComputed = n;
        }
    }

    if (lastComputed == GenTraits<gen>::power2) {
        std::cout << "initializing from base matrix: F^(2^" << GenTraits<gen>::power2 << ")\n";
    }
    else {
        std::string fn = mkFileName(filepath, lastComputed);
        std::cout << "initializing F^(2 ^ " << lastComputed << ") from file " << fn << "\n";
        std::ifstream is(fn, ios::binary);
        f[lastComputed % 2].fromBin(is);
    }

    std::cout << "  ";
    f[lastComputed % 2].printSparsity();

    std::vector<typename GenTraits<gen>::matrix_t::buffer_t> buffers(nThreads);

    for (size_t i = lastComputed + 1; i <= 19937 && (lastComputed < stopIndex); ++i) {
        std::cout << "computing F^(2^" << i << ")\n";
        size_t in = (i + 1) % 2;
        size_t out = (i) % 2;
        f[out].resetZero();
        square(f[in], f[out], buffers);
        std::cout << "  ";
        f[out].printSparsity();

        if ((i % saveFrequency) == 0 || i > 19930) {
            std::string fn = mkFileName(filepath, i);
            std::cout << "  saving file: " << fn << " ... ";
            std::ofstream of(fn, ios::binary);
            f[out].toBin(of);
            of.close();
            std::cout << "saved\n";
            lastComputed = i;
        }

    }
}


int main(int argc, const char** argv)
{
    // parse command line arguments
    size_t nThreads = std::thread::hardware_concurrency();
    std::string filepath = "./dat/";
    std::string gentype;
    size_t saveFrequency = 100;
    size_t stopIndex = std::numeric_limits<size_t>::max();
    if (argc % 2 == 0)
        usage();
    for (int i = 1; i < argc; i += 2) {
        string key(argv[i]);
        const char* value = argv[i + 1];
        if (key == "-j") {
            nThreads = atoi(value);
        }
        else if (key == "-p") {
            filepath = value;
            if (filepath.back() != '/')
                filepath.push_back('/');
        }
        else if (key == "-g") {
            gentype = value;
        }
        else if (key == "-f") {
            saveFrequency = atoi(value);
        }
        else if (key == "-s") {
            stopIndex = atoi(value);
        }
        else
            usage();
    }

    std::cout << "nThreads = " << nThreads << "\n";
    std::cout << "filepath = " << filepath << "\n";
    std::cout << "savefrequency = " << saveFrequency << "\n";
    std::cout << "generator = " << gentype << "\n";
    std::cout << "stopindex = " << stopIndex << "\n";

    if (gentype == "mt")
        run<mt>(filepath, nThreads, saveFrequency, stopIndex);
    if (gentype == "sfmt")
        run<sfmt>(filepath, nThreads, saveFrequency, stopIndex);
    else
        usage();

    return 0;

}
