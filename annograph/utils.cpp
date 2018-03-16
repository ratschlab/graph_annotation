#include "utils.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <cassert>

namespace utils {

/**
 *  Returns the input file type, given a filename
 */
std::string get_filetype(const std::string &fname) {
    size_t dotind = fname.rfind(".");
    if (dotind == std::string::npos)
        return "";

    std::string ext = fname.substr(dotind);

    if (ext == ".gz") {
        size_t nextind = fname.substr(0, dotind - 1).rfind(".");
        if (nextind == std::string::npos)
            return "";

        ext = fname.substr(nextind, dotind - nextind);
    }

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".vcf") {
        return "VCF";
    } else if ((ext == ".fq") || (ext == ".fastq")) {
        return "FASTQ";
    } else {
        return "FASTA";
    }
}

/**
 * Given a minimum number of splits,
 * generate a list of suffices from the alphabet.
 */
std::deque<std::string> generate_strings(const std::string &alphabet,
                                         size_t length) {

    std::deque<std::string> suffices = { "" };
    while (suffices[0].length() < length) {
        for (const char c : alphabet) {
            suffices.push_back(c + suffices[0]);
        }
        suffices.pop_front();
    }
    assert(suffices.size() == std::pow(alphabet.size(), length));
    return suffices;
}

//These next functions borrowed from libmaus2
uint64_t deserializeNumber(std::istream &in) {
   uint64_t n = 0;
    for (size_t i = 0; i < sizeof(n); ++i) {
        n = (n << 8) | in.get();
    }
    return n;
}

uint64_t serializeNumber(std::ostream &out, uint64_t const n) {
    out.put((n >> (7 * 8)) & 0xFF);
    out.put((n >> (6 * 8)) & 0xFF);
    out.put((n >> (5 * 8)) & 0xFF);
    out.put((n >> (4 * 8)) & 0xFF);
    out.put((n >> (3 * 8)) & 0xFF);
    out.put((n >> (2 * 8)) & 0xFF);
    out.put((n >> (1 * 8)) & 0xFF);
    out.put((n >> (0 * 8)) & 0xFF);

    if (!out) {
        std::cerr << "Serialization failure" << std::endl;
        exit(1);
    }

    return 8;
}

uint64_t serializeNumberVector(std::ostream &out, std::vector<uint64_t> const &v) {
    uint64_t  s = 0;
    s += serializeNumber(out, v.size());
    for (auto &n : v) {
        s += serializeNumber(out, n);
    }
    return s;
}

} // namespace utils
