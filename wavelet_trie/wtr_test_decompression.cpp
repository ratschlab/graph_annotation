#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <omp.h>
#include <thread>
#include <future>
#include <mutex>
#include <chrono>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>

#include "wavelet_trie.hpp"
#include "unix_tools.hpp"
#include "getRSS.h"


//TODO: replace getenv with real argument parsing

typedef annotate::cpp_int cpp_int;


int main(int argc, char** argv) {
    if (argc <= 1) {
        std::cerr << "ERROR: please pass in at least one file." << std::endl;
        exit(1);
    }
    const char *maxcount = std::getenv("MAX_COUNT");
    size_t max_count = 0;
    if (maxcount) {
        max_count = atol(maxcount);
    }

    //environment variable arguments

    std::vector<std::future<annotate::WaveletTrie*>> wtrs;
    wtrs.reserve(argc);

    std::cout << "Starting merge\n";

    for (int i = 1; i < argc; ++i) {
        std::cout << "Pushing graph " << i - 1 << "\n";
        wtrs.push_back(std::async(std::launch::deferred, [i, &argv]() {
            std::ifstream fin(argv[i]);
            annotate::WaveletTrie *wtr = new annotate::WaveletTrie();
            wtr->load(fin);
            fin.close();
            return wtr;
        }));
    }

    std::cout << "Loading graph\n";
    auto *wtr = wtrs[0].get();

    if (wtr) {
        if (!max_count) {
            max_count = wtr->size();
        } else {
            max_count = std::min(max_count, wtr->size());
        }
        Timer extract_timer;
        std::cout << "Input:" << std::endl;
        std::cout << "Num edges:\t" << wtr->size() << std::endl;
        extract_timer.reset();
        for (size_t i = 0; i < max_count; ++i) {
            wtr->at(i);
        }
        double query_time = extract_timer.elapsed();
        std::cout << "Decomptime:\t" << query_time << std::endl;
        std::cout << "Query time:\t" << query_time / wtr->size() << std::endl;
        delete wtr;
    }
    return 0;
}


