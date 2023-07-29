#include "jump_matrix.h"

#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;

enum class EncodeMode {
    Base64,
    Binary,
    Hex,
    Header
};

auto parseFilePath(const std::string& filePath)
{
    std::string folder, filename, extension;
    EncodeMode mode = EncodeMode::Header; // Default mode

    // Find the last occurrence of the path separator '/' to separate folder and filename.
    size_t lastSeparatorPos = filePath.find_last_of('/');
    if (lastSeparatorPos != std::string::npos) {
        folder = filePath.substr(0, lastSeparatorPos); // Extract the folder part.
        filename = filePath.substr(lastSeparatorPos + 1); // Extract the filename part.
    }
    else {
        // If there's no separator, the file is in the current folder or it has no folder part.
        folder = "";
        filename = filePath; // The whole filePath is the filename in this case.
    }

    // Find the last occurrence of the dot '.' to separate filename and extension.
    size_t lastDotPos = filename.find_last_of('.');
    if (lastDotPos != std::string::npos) {
        extension = filename.substr(lastDotPos + 1); // Extract the extension part.
        filename = filename.substr(0, lastDotPos); // Update the filename without the extension.
    }

    // Validate the extension. Set the EncodeMode accordingly.
    if (extension == "b64") {
        mode = EncodeMode::Base64;
    }
    else if (extension == "bits") {
        mode = EncodeMode::Binary;
    }
    else if (extension == "hex") {
        mode = EncodeMode::Hex;
    }
    else if (extension == "h") {
        mode = EncodeMode::Header;
    }
    else {
        THROW("Invalid file extension in filename " << filePath << ". Extension should be one of 'b64', 'bits', 'hex', or 'h'.");
    }

    return std::make_tuple(folder, filename, extension, mode);
}
/*
template <size_t nRows, size_t nCols>
void testEncoder(const BinaryMatrix<nRows, nCols>& m, EncodeMode enc)
{
    const char* modename = enc == Base64 ? "base64" : "hex";

    BinaryMatrix<nRows, nCols> m2;

    std::cout << "saving matrix to " << modename << " stream\n";
    std::ostringstream os;
    if (enc == Base64)
        m.toBase64(os);
    else
        m.toHex(os);

    std::cout << "first 16 characters of the stream\n";
    std::string s = os.str();
    for (size_t i = 0; i < 16; ++i)
        std::cout << s[i];
    std::cout << "\n";

    std::cout << "reading back the matrix from " << modename << " stream\n";
    std::istringstream is(os.str());
    if (enc == Base64)
        m2.fromBase64(is);
    else
        m2.fromHex(is);

    std::cout << "compare with original matrix\n";
    if (!(m == m2))
        throw std::invalid_argument("error in roundtrip");

    std::cout << "completed\n";
}
*/
void usage()
{
    std::cerr
        << "Invalid command line arguments\n"
        << "Syntax:\n"
        << "   encoder -i inputfile -o outputfile\n"
        << "Example:\n"
        << "   encoder -i dir1/F19937.bits -o dir2/F19937.b64\n"
        << "valid file extensins are: b64, hex, bits, h\n";
    THROW("");
}


int main(int argc, const char** argv)
{
    // parse command line arguments
    string inputfile, outputfile;
    if (argc % 2 == 0)
        usage();
    for (int i = 1; i < argc; i += 2) {
        string key(argv[i]);
        string value(argv[i + 1]);
        if (key == "-i")
            inputfile = value;
        else if (key == "-o")
            outputfile = value;
        else {
            usage();
        }
    }

    // parse filenames
    auto [ifolder, iname, iext, imode] = parseFilePath(inputfile);
    auto [ofolder, oname, oext, omode] = parseFilePath(outputfile);

    std::cout << "input: " << inputfile << ", output: " << outputfile << ": ... ";

    try {
        MT19937Matrix matrix;

        std::ifstream is;
        if (imode == EncodeMode::Hex) {
            is.open(inputfile);
            MYASSERT(is.is_open(), "error opening binary file: " << inputfile);
            matrix.fromHex(is);
        }
        if (imode == EncodeMode::Base64) {
            is.open(inputfile);
            MYASSERT(is.is_open(), "error opening binary file: " << inputfile);
            matrix.fromBase64(is);
        }
        else if (imode == EncodeMode::Binary) {
            is.open(inputfile, std::ios::binary);
            MYASSERT(is.is_open(), "error opening binary file: " << inputfile);
            matrix.fromBin(is);
        }
        else {
            THROW("input file type not supported");
        }

        std::ofstream os;
        if (omode == EncodeMode::Base64) {
            os.open(outputfile);
            MYASSERT(os.is_open(), "error opening binary file: " << outputfile);
            matrix.toBase64(os);
        }
        else if (omode == EncodeMode::Hex) {
            os.open(outputfile);
            MYASSERT(os.is_open(), "error opening binary file: " << outputfile);
            matrix.toHex(os);
        }
        else if (omode == EncodeMode::Header) {
            os.open(outputfile);
            MYASSERT(os.is_open(), "error opening binary file: " << outputfile);
            matrix.toArrayChar(os);
        }
        else if (omode == EncodeMode::Binary) {
            os.open(outputfile, std::ios::binary);
            MYASSERT(os.is_open(), "error opening binary file: " << outputfile);
            matrix.toBin(os);
        }
        else {
            THROW("output file type not supported");
        }
        MYASSERT(os.is_open(), "error opening binary file: " << outputfile);


    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return -1;
    }
    catch (...) {
        std::cerr << "Error: unknown exception\n";
        return -1;
    }
    std::cout << "done\n";
    return 0;
}