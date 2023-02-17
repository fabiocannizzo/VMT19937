#include "jump_matrix.h"

#include <thread>
#include <memory>
#include <sstream>

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

int main(int narg, const char *args[])
{
    MT19937Matrix f[2];
    f[0].init();
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


    for (size_t i = 0; i < 19936; ++i) {
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

        std::ofstream os("temp");
        f[out].toBase64(os);
        os.close();
        std::ifstream is("temp");
        f[in].resetZero();
        f[in].fromBase64(is);
        is.close();
        if (!(f[in] == f[out]))
            std::cout << "read save problem\n";
/*
        if ((i % 200) == 0 || i > 19930) {
            std::ostringstream os;
            os << "F" << std::setw(5) << std::setfill('0') << i+1;
            f[out].base64Out("", os.str().c_str());
        }
*/

    }

    return 0;
}
