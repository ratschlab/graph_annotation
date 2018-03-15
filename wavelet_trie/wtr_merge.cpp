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
size_t total_bits = 0;
size_t set_bits = 0;


int main(int argc, char** argv) {
    std::cout << "sizeof(Node):\t" << sizeof(annotate::WaveletTrie::Node) << "\n";
    std::cout << "sizeof(alpha_t):\t" << sizeof(annotate::alpha_t) << "\n";
    std::cout << "sizeof(beta_t):\t" << sizeof(annotate::beta_t) << "\n";
    std::cout << "sizeof(rrr_t):\t" << sizeof(annotate::rrr_t) << "\n";
    std::cout << "sizeof(sd_t):\t" << sizeof(annotate::sd_t) << "\n";
    std::cout << "sizeof(bv_t:rank1):\t" << sizeof(annotate::rank1_t) << "\n";
    std::cout << "sizeof(bv_t:rank0):\t" << sizeof(annotate::rank0_t) << "\n";
    std::cout << "sizeof(children):\t" << sizeof(annotate::WaveletTrie::Node*) << "\n";
    std::cout << "sizeof(cpp_int):\t" << sizeof(cpp_int) << "\n";
    std::cout << "sizeof(mpz_t):\t" << sizeof(mpz_t) << "\n";
    std::cout << "sizeof(bv_t):\t" << sizeof(annotate::bv_t) << "\n";
    if (argc <= 3) {
        std::cerr << "ERROR: please pass in at least three files: two inputs and an output." << std::endl;
        exit(1);
    }

    //environment variable arguments
    const char *njobs = std::getenv("NJOBS");
    size_t n_jobs = 1;
    if (njobs) {
        n_jobs = atoi(njobs);
        omp_set_num_threads(n_jobs);
    }

    std::vector<std::future<annotate::WaveletTrie*>> wtrs;
    wtrs.reserve(argc - 2);

    std::cout << "Starting merge\n";

    for (int i = 1; i < argc - 1; ++i) {
        std::cout << "Pushing graph " << i - 1 << "\n";
        wtrs.push_back(std::async(std::launch::deferred, [i, &argv]() {
            std::ifstream fin(argv[i]);
            annotate::WaveletTrie *wtr = new annotate::WaveletTrie();
            wtr->load(fin);
            fin.close();
            return wtr;
        }));
    }
/*
    std::cout << "Loading graph 0\n";
    wtrs[0].wait_for(std::chrono::seconds(1));
    annotate::WaveletTrie *curwtr = NULL;
    if (wtrs.size() > 1) {
        std::cout << "Loading graph 1\n";
        wtrs[1].wait_for(std::chrono::seconds(1));
    }

    auto *wtr = wtrs[0].get();
    std::cout << "Graph 0 loaded\n";
    curwtr = wtrs[1].get();
    std::cout << "Graph 1 loaded\n";
    for (auto it = wtrs.begin() + 2; it != wtrs.end(); ++it) {
        std::cout << "Loading graph " << std::distance(wtrs.begin(), it) << "\n";
        it->wait_for(std::chrono::seconds(1));
        std::cout << "Merging graph " << (std::distance(wtrs.begin(), it) - 1) << "\n";
        wtr->insert(*curwtr);
        delete curwtr;
        std::cout << "." << std::flush;
        curwtr = it->get();
        std::cout << "Graph " << std::distance(wtrs.begin(), it) << " loaded\n";
    }
    if (curwtr) {
        std::cout << "Merging graph " << (std::distance(wtrs.begin(), wtrs.end()) - 1) << "\n";
        wtr->insert(*curwtr);
        delete curwtr;
        std::cout << "." << std::flush;
    }
    */
    std::cout << "Loading graph 0\n";
    auto *wtr = wtrs[0].get();
    for (auto it = wtrs.begin() + 1; it != wtrs.end(); ++it) {
        std::cout << "Loading graph " << std::distance(wtrs.begin(), it) << "\n";
        auto *curwtr = it->get();
        wtr->insert(*curwtr);
        delete curwtr;
    }
    std::cout << std::endl;

    if (wtr) {
        std::cout << "Serializing:\t" << std::flush;
        std::ofstream fout(argv[argc - 1], std::ofstream::binary | std::ofstream::ate);
        if (!fout.good()) {
            std::cerr << "ERROR: bad file " << argv[argc - 1] << std::endl;
            exit(1);
        }
        auto stats = wtr->serialize(fout);
        std::cout << "Input:" << std::endl;
        std::cout << "Num edges:\t" << wtr->size() << std::endl;
        std::cout << "Total bits:\t" << total_bits << std::endl;
        std::cout << "Set bits:\t" << set_bits << std::endl;
        std::cout << "Wavelet trie:" << std::endl;
        std::cout << "Num leaves:\t" << stats << std::endl;
        std::cout << "Num bytes:\t" << fout.tellp() << std::endl;
        std::cout << "Bits/edge:\t" << (double)fout.tellp() * 8.0 / (double)wtr->size() << std::endl;
        fout.close();
        delete wtr;
    }
    std::cout << "Done\n";
    return 0;
}


