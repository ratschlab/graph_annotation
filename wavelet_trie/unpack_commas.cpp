#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>
#include "wavelet_trie.hpp"

int main(int argc, char **argv) {
    if (argc < 3)
        exit(1);
    std::ifstream fin(argv[1]);
    std::ofstream fout(argv[2]);

    std::unordered_map<std::string, std::set<size_t>> inds;
    boost::archive::binary_iarchive ia(fin);
    ia & inds;
    for (auto &kmer : inds) {
        auto it = kmer.second.begin();
        if (kmer.second.size()) {
            fout << *it;
            ++it;
        }
        for (; it != kmer.second.end(); ++it) {
            fout << "," << *it;
        }
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}
