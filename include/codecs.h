#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cctype>
#include <exception>

using namespace std;

namespace Encoder {

    static const string base64Chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    static const string hexChars =
        "0123456789"
        "ABCDEF";

    typedef unsigned char uint8_t;

    inline bool isBase64(uint8_t c) {
        return (isalnum(c) || (c == '+') || (c == '/'));
    }

    inline bool isHex(uint8_t c) {
        return (isdigit(c) || ((c >= 'a') && (c <= 'f')) || ((c >= 'A') && (c <= 'F')));
    }


    inline uint8_t hexToDec(char hex)
    {
        if (hex >= '0' && hex <= '9')
            return hex - '0';
        if (hex >= 'a' && hex <= 'f')
            return hex - 'a' + 10;
        if (hex >= 'A' && hex <= 'F')
            return hex - 'A' + 10;
        char t[2] = { hex, 0 };
        cout << "not a hex character: " << t << endl;
        throw std::invalid_argument("not a hex character");
    }

    inline uint8_t b64ToDec(char b64)
    {
        const size_t nAlpha = 'Z' - 'A' + 1;
        if (b64 >= 'A') {
            if (b64 <= 'Z')
                return b64 - 'A';
            if (b64 >= 'a' && b64 <= 'z')
                return b64 - 'a' + nAlpha;
        }
        if (b64 >= '0' && b64 <= '9')
            return b64 - '0' + 2 * nAlpha;
        if (b64 == '+')
            return 62;
        if (b64 == '/')
            return 63;
        char t[2] = { b64, 0 };
        cout << "not a base 64 character: " << t << endl;
        throw std::invalid_argument("not a b64 character");
    }

    inline uint8_t hexPairToDec(const char hex[2])
    {
        return (hexToDec(hex[0]) << 4) | hexToDec(hex[1]);
    }


    inline void hexToText(string& res, const string& hex)
    {
        size_t n = hex.size();
        size_t m = n / 2;
        assert((n & 0x1) == 0);  // n must be even
        res.resize(m);
        for (size_t i = 0, j = 0; i < m; ++i, j += 2)
            res[i] = hexPairToDec(&hex[j]);
    }

    inline void textToHex(string& hex, const string& text)
    {
        size_t n = text.size();
        size_t m = n * 2;
        const uint8_t* p = (const uint8_t*) text.c_str();
        hex.resize(m);
        for (size_t i = 0, j = 0; i < n; ++i, j += 2) {
            hex[j] = hexChars[p[i] >> 4];
            hex[j + 1] = hexChars[p[i] & 0xF];
        }
    }


    inline void dec3ToBase64(char res[4], const uint8_t dec[3])
    {
        res[0] = base64Chars[dec[0] >> 2];
        res[1] = base64Chars[((dec[0] & 0x3) << 4) | (dec[1] >> 4)];
        res[2] = base64Chars[((dec[1] & 0xF) << 2) | (dec[2] >> 6)];
        res[3] = base64Chars[dec[2] & 0x3F];
    }

    inline void base64ToDec3(uint8_t res[3], const uint8_t b64[4])
    {
        const uint8_t dec[4] = { b64ToDec(b64[0]), b64ToDec(b64[1]), b64ToDec(b64[2]), b64ToDec(b64[3]) };
        res[0] = (dec[0] << 2) | (dec[1] >> 4);
        res[1] = (dec[1] << 4) | (dec[2] >> 2);
        res[2] = (dec[2] << 6) | dec[3];
    }


    inline void textToBase64(string& res, const string& txt)
    {
        size_t n = txt.size();
        res.resize((n / 3 + ((n % 3) ? 1 : 0)) * 4);
        const uint8_t* dec = reinterpret_cast<const uint8_t*>(txt.c_str());
        size_t j = 0, i = 0;
        for (; i + 2 < n; i += 3, j += 4)
            dec3ToBase64(&res[j], dec + i);
        switch (n % 3) {
            case 0:
                break;
            case 1:
            {
                uint8_t aux[3] = { dec[i], 0, 0 };
                dec3ToBase64(&res[j], aux);
                res[j + 2] = res[j + 3] = '=';
            }
            break;
            case 2:
            {
                uint8_t aux[3] = { dec[i], dec[i + 1], 0 };
                dec3ToBase64(&res[j], aux);
                res[j + 3] = '=';
            }
            break;
        };
    }

    inline void base64ToText(string& res, const string& b64)
    {
        size_t n = b64.size();
        if (!((n % 4) == 0 && n > 0)) {
            std::cout << "bad stream size: " << n << "\n";
            throw std::invalid_argument("bad stream size");
        }
        //assert((n % 4) == 0 && n > 0);
        size_t nEq = (b64[n - 1] == '=' ? 1 : 0) + (b64[n - 2] == '=' ? 1 : 0);
        size_t n4 = n / 4 - (nEq > 0);
        res.resize(n4 * 3 + ((nEq > 0) ? 3 - nEq : 0));
        const uint8_t* dec = reinterpret_cast<const uint8_t*>(b64.c_str());
        size_t i = 0, j = 0;
        for (; i < n4; ++i, j += 3)
            base64ToDec3(reinterpret_cast<uint8_t*>(&res[j]), dec + 4 * i);
        switch (nEq) {
            case 0:
                break;
            case 1:
            {
                uint8_t r[3];
                const uint8_t b[4] = { dec[4 * i],dec[4 * i + 1],dec[4 * i + 2],'A' };
                base64ToDec3(r, b);
                res[j++] = r[0];
                res[j] = r[1];
            }
            break;
            case 2:
            {
                uint8_t r[3];
                const uint8_t b[4] = { dec[4 * i],dec[4 * i + 1],'A','A' };
                base64ToDec3(r, b);
                res[j] = r[0];
            }
            break;
        };
    }

};
