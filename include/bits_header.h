#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

namespace xvmt {

enum class BitsFileType : uint8_t {
    Polynomial = 0,
    Matrix = 1
};

enum class BitsGenType : uint8_t {
    MT32 = 0,
    MT64 = 1,
    SFMT = 2,
    Unknown = 255
};

inline const char* toString(BitsGenType g) {
    switch (g) {
        case BitsGenType::MT32: return "mt32";
        case BitsGenType::MT64: return "mt64";
        case BitsGenType::SFMT: return "sfmt";
        default: return "unknown";
    }
}

inline BitsGenType stringToGenType(const std::string& s) {
    if (s == "mt32") return BitsGenType::MT32;
    if (s == "mt64") return BitsGenType::MT64;
    if (s == "sfmt") return BitsGenType::SFMT;
    return BitsGenType::Unknown;
}

#pragma pack(push, 1)
struct BitsHeader {
    static constexpr char MAGIC_STR[] = "VMTBITS";
    static constexpr uint8_t VERSION_NUM = 1;

    char magic[8];         // "VMTBITS\0"
    uint8_t version;       // 1
    uint8_t fileType;      // BitsFileType
    uint8_t genType;       // BitsGenType
    uint32_t jumpPower2;   // N where jump is 2^N
    uint64_t rows;
    uint64_t cols;
    uint8_t reserved[33];  // Pad to 64 bytes total

    BitsHeader() {
        std::memset(this, 0, sizeof(BitsHeader));
        std::memcpy(magic, MAGIC_STR, 7); // "VMTBITS"
        version = VERSION_NUM;
    }

    bool isValid() const {
        return std::memcmp(magic, MAGIC_STR, 7) == 0 && version == VERSION_NUM;
    }

    void write(std::ostream& os) const {
        os.write(reinterpret_cast<const char*>(this), sizeof(BitsHeader));
    }

    bool read(std::istream& is) {
        auto pos = is.tellg();
        is.read(reinterpret_cast<char*>(this), sizeof(BitsHeader));
        if (is.gcount() == sizeof(BitsHeader) && isValid()) {
            return true;
        }
        is.clear();
        is.seekg(pos);
        return false;
    }
};
#pragma pack(pop)

static_assert(sizeof(BitsHeader) == 64, "BitsHeader size must be 64 bytes");

} // namespace xvmt
