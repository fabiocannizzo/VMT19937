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

typedef BinarySquareMatrix<19937>  MT19937Matrix;

// initialize the matrix as per MT19937 32 bit generartor transition matrix
void initMT19937(MT19937Matrix& m)
{
    static_assert(MT19937Matrix::s_nBits == 19937);
    static const size_t s_nBits = MT19937Matrix::s_nBits;
    static const size_t s_nWordBits = MT19937Matrix::s_nWordBits;
    static const uint32_t s_matA = 0x9908B0DF;
    static const uint32_t s_M = 397;

    // from row 0 to to row nBits - 32, state bits are just shifted left by 32 bits
    for (uint32_t r = 0; r < s_nBits - s_nWordBits; ++r)
        m.setBit(r, r + s_nWordBits);

    // the new state element is composed of element which was in position M
    for (uint32_t i = 0; i < s_nWordBits; i++)
        m.setBit(s_nBits - s_nWordBits + i, 1 + (s_M - 1) * s_nWordBits + i);
    for (uint32_t i = 0; i < s_nWordBits; ++i)
        if (s_matA & (uint32_t(1) << i))
            m.setBit(s_nBits - s_nWordBits + i, 1);
    m.setBit(s_nBits - 2, 0);
    for (uint32_t i = 0; i < s_nWordBits - 2; ++i)
        m.setBit(s_nBits - s_nWordBits + i, 2 + i);
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

std::string mkFileName(size_t n)
{
    std::ostringstream os;
    os << "F" << std::setw(5) << std::setfill('0') << n << extension;
    return os.str();
}

int main(int narg, const char* args[])
{
    MT19937Matrix f[2];

    int lastComputed = -1;

    for (const auto& entry : std::filesystem::directory_iterator("./")) {
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
        std::string fn = mkFileName(lastComputed);
        std::cout << "initializing from file "<< fn << "\n";
        std::ifstream is(fn);
        f[lastComputed % 2].fromBase64(is);
    }

    std::cout << "  ";
    f[lastComputed % 2].printSparsity();

    const size_t nThreads = 4;
    const size_t saveFrequency = 100;

    std::vector<MT19937Matrix::buffer_t> buffers(nThreads);

    for (size_t i = lastComputed + 1; i <= 19937; ++i) {
        std::cout << "computing F^(2^" << i << ")\n";
        size_t in = (i + 1) % 2;
        size_t out = (i) % 2;
        f[out].resetZero();
        square(f[in], f[out], buffers, nullptr);
        std::cout << "  ";
        f[out].printSparsity();

        if ((i  % saveFrequency) == 0 || i > 19930) {
            std::string fn = mkFileName(i);
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
