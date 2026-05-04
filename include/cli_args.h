#pragma once

#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <type_traits>

using ArgMap = std::map<std::string, std::string>;

inline ArgMap parseArgs(int argc, const char** argv)
{
    ArgMap args;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.empty()) continue;
        size_t start = 0;
        if (s.size() >= 2 && s[0] == '-' && s[1] == '-') start = 2;
        else if (s.size() >= 1 && s[0] == '-') start = 1;
        else continue;

        std::string body = s.substr(start);
        size_t eq = body.find('=');
        if (eq != std::string::npos) {
            args[body.substr(0, eq)] = body.substr(eq + 1);
        }
        else {
            args[body] = "";
        }
    }
    return args;
}

template <typename T>
inline bool consumeArg(ArgMap& args, const std::string& key, bool compulsory, T& result)
{
    auto it = args.find(key);
    if (it != args.end()) {
        std::string val = it->second;
        args.erase(it);

        if constexpr (std::is_same_v<T, std::string>) {
            result = val;
        }
        else {
            if (val.empty()) {
                throw std::runtime_error("Argument " + key + " requires a value");
            }
            std::stringstream ss(val);
            if (!(ss >> result) || !ss.eof()) {
                throw std::runtime_error("Argument conversion failed for key: " + key + " with value: " + val);
            }
        }
        return true;
    }

    if (compulsory) {
        throw std::runtime_error("Compulsory argument missing: " + key);
    }

    return false;
}

inline bool consumeArg(ArgMap& args, const std::string& key)
{
    auto it = args.find(key);
    if (it != args.end()) {
        args.erase(it);
        return true;
    }
    return false;
}

inline void waitForDebugger(ArgMap& args)
{
    if (consumeArg(args, "w") || consumeArg(args, "wait")) {
        std::cout << "press a key to continue..." << std::endl;
        std::cin.get();
    }
}
