#include "polynomial_jump.h"
#include "jump_matrix.h"
#include "Params.h"
#include <filesystem>
#include <iostream>
#include <string>
#include <regex>
#include <fstream>

using namespace std;
using namespace xvmt;
using namespace xvmt::details;
namespace fs = std::filesystem;

void migrate_poly(const fs::path& path, BitsGenType genType) {
    // Check if it already has a header
    {
        ifstream ifs(path, ios::binary);
        BitsHeader h;
        if (h.read(ifs)) {
            cout << "Skipping already migrated poly: " << path.string() << endl;
            return;
        }
    }

    cout << "Migrating poly: " << path.string() << endl;
    Polynomial<32768> p;
    {
        ifstream ifs(path, ios::binary);
        p.fromBinStream(ifs);
    }
    
    // Infer jump from filename J<N>
    uint32_t jump = 0;
    string name = path.filename().string();
    regex re("J(\\d+)");
    smatch match;
    if (regex_search(name, match, re)) {
        jump = (uint32_t)stoi(match[1]);
    }
    
    {
        ofstream ofs(path, ios::binary);
        p.toBin(ofs, genType, jump);
    }
}

template<typename Matrix>
void migrate_matrix(const fs::path& path, BitsGenType genType) {
    // Check if it already has a header
    {
        ifstream ifs(path, ios::binary);
        BitsHeader h;
        if (h.read(ifs)) {
            cout << "Skipping already migrated matrix: " << path.string() << endl;
            return;
        }
    }

    cout << "Migrating matrix: " << path.string() << endl;
    Matrix m;
    {
        ifstream ifs(path, ios::binary);
        m.fromBinStream(ifs);
    }
    
    // Infer jump from filename F<N>
    uint32_t jump = 0;
    string name = path.filename().string();
    regex re("F(\\d+)");
    smatch match;
    if (regex_search(name, match, re)) {
        jump = (uint32_t)stoi(match[1]);
    }
    
    {
        ofstream ofs(path, ios::binary);
        m.toBin(ofs, genType, jump);
    }
}

int main() {
    try {
        // Poly mt32
        if (fs::exists("dat/poly/mt32")) {
            for (const auto& entry : fs::directory_iterator("dat/poly/mt32")) {
                if (entry.path().extension() == ".bits") migrate_poly(entry.path(), BitsGenType::MT32);
            }
        }
        // Poly mt64
        if (fs::exists("dat/poly/mt64")) {
            for (const auto& entry : fs::directory_iterator("dat/poly/mt64")) {
                if (entry.path().extension() == ".bits") migrate_poly(entry.path(), BitsGenType::MT64);
            }
        }
        // Poly sfmt
        if (fs::exists("dat/poly/sfmt")) {
            for (const auto& entry : fs::directory_iterator("dat/poly/sfmt")) {
                if (entry.path().extension() == ".bits") migrate_poly(entry.path(), BitsGenType::SFMT);
            }
        }
        
        // Matrix mt32
        if (fs::exists("dat/matrix/mt32")) {
            for (const auto& entry : fs::directory_iterator("dat/matrix/mt32")) {
                if (entry.path().extension() == ".bits") migrate_matrix<MT19937Matrix<32>>(entry.path(), BitsGenType::MT32);
            }
        }
        // Matrix mt64
        if (fs::exists("dat/matrix/mt64")) {
            for (const auto& entry : fs::directory_iterator("dat/matrix/mt64")) {
                if (entry.path().extension() == ".bits") migrate_matrix<MT19937Matrix<64>>(entry.path(), BitsGenType::MT64);
            }
        }
        // Matrix sfmt
        if (fs::exists("dat/matrix/sfmt")) {
            for (const auto& entry : fs::directory_iterator("dat/matrix/sfmt")) {
                if (entry.path().extension() == ".bits") migrate_matrix<SFMT19937Matrix>(entry.path(), BitsGenType::SFMT);
            }
        }
    } catch (const exception& e) {
        cerr << "Error during migration: " << e.what() << endl;
        return 1;
    }
    
    cout << "Migration complete." << endl;
    return 0;
}
