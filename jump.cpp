#include "jump_matrix.h"

#include <thread>
#include <memory>
#include <sstream>
#include <filesystem>

void wait()
{
    std::cout << "press a key...";
    char c;
    std::cin >> c;
}

void square(const MT19937Matrix& src, MT19937Matrix& dst, std::vector<MT19937Matrix::buffer_t>& buffers, const MT19937Matrix* target)
{
    auto start = std::chrono::system_clock::now();
    dst.square(src, buffers, target);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "done in: " << std::fixed << std::setprecision(2) << elapsed_seconds.count() << "s" << std::endl;
}

const std::string extension = ".b64";

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
        << "jump [-j nthreads] [-p filepath] [-f savefrequency]\n"
        << "  nthreads defaults to 1\n"
        << "  filepath defaults to ./"
        << "  savefrequency defaults to 100";
    std::exit(-1);
}

int main(int argc, const char** argv)
{
    // parse command line arguments
    size_t nThreads = std::thread::hardware_concurrency();
    std::string filepath = "./dat/";
    size_t saveFrequency = 100;
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
        else if (key == "-f") {
            saveFrequency = atoi(value);
        }
        else
            usage();
    }
    std::cout << "nThreads = " << nThreads << "\n";
    std::cout << "filepath = " << filepath << "\n";
    std::cout << "savefrequency = " << saveFrequency << "\n";

    MT19937Matrix f[2];

    int lastComputed = -1;

    for (const auto& entry : std::filesystem::directory_iterator(filepath)) {
        if (std::filesystem::is_regular_file(entry) && entry.path().has_extension() && entry.path().extension().string() == extension) {
            std::string s = entry.path().filename().string();
            s = s.substr(1, s.length() - 1 - extension.length());
            int n = atoi(s.c_str());
            if (n > lastComputed)
                lastComputed = n;
        }
    }

    if (lastComputed == -1) {
        std::cout << "initializing from base matrix\n";
        initMT19937(f[0]);
        lastComputed = 0;
    }
    else {
        std::string fn = mkFileName(filepath, lastComputed);
        std::cout << "initializing from file " << fn << "\n";
        std::ifstream is(fn);
        f[lastComputed % 2].fromBase64(is);
    }

    std::cout << "  ";
    f[lastComputed % 2].printSparsity();

    std::vector<MT19937Matrix::buffer_t> buffers(nThreads);

    for (size_t i = (size_t) lastComputed + 1; i <= 19937; ++i) {
        std::cout << "computing F^(2^" << i << ")\n";
        size_t in = (i + 1) % 2;
        size_t out = (i) % 2;
        f[out].resetZero();
        square(f[in], f[out], buffers, nullptr);
        std::cout << "  ";
        f[out].printSparsity();

        if ((i  % saveFrequency) == 0 || i > 19930) {
            std::string fn = mkFileName(filepath, i);
            std::cout << "  saving file: " << fn << " ... ";
            std::ofstream of(fn);
            f[out].toBase64(of);
            std::cout << "saved\n";
        }

    }

    return 0;

}

#if 0
int main(int narg, const char *args[])
{

    MT19937Matrix f[2];
    initMT19937(f[0]);
    f[0].printBits(MT19937Matrix::s_nBits - 32, 0, 32, 32);
    std::cout << "f0: ";  f[0].printSparsity();

    MT19937Matrix m[3];
    for (size_t i = 0; i < 3; ++i) {
        char n[] = { 'm', (char) (i + '0'), 0 };
        std::ostringstream os;
        os << "r" << i << ".txt";
        m[i].fromMatlabSparseFile(os.str().c_str());
        std::cout << n << ": ";
        m[i].printSparsity();
    }

    if (!(f[0] == m[0]))
        std::cout << "base matrix is different\n";

    const size_t nThreads = 4;

    std::vector<MT19937Matrix::buffer_t> buffers(nThreads);


    for (size_t i = 0; i <= 19937; ++i) {
        std::cout << "computing F^(2^" << i + 1 << ")\n";
        size_t in = i % 2;
        size_t out = (i + 1) % 2;
        f[out].resetZero();
//        if (i == 0)
            square(f[in], f[out], buffers, nullptr);
//        else
//            square(f[in], f[out], buffers, &target);
        //if (i==0)
        //    f[1].matlabOut("xR2");
        if (i < 2)
            if (!(f[out] == m[i+1]))
                std::cout << "base matrix is different\n";
        f[out].printSparsity();

#if 0
        std::ofstream os("temp");
        f[out].toBase64(os);
        os.close();
        std::ifstream is("temp");
        f[in].resetZero();
        f[in].fromBase64(is);
        is.close();
        if (!(f[in] == f[out]))
            std::cout << "read save problem\n";
#endif

        if ((i % 200) == 0 || i > 19930) {
            std::ostringstream os;
            os << "F" << std::setw(5) << std::setfill('0') << i+1;
            std::cout << "saving file: " << os.str() << "...";
            std::ofstream of(os.str());
            f[out].toBase64(of);
            std::cout << "saved\n";
        }


    }

    return 0;
}

#endif
