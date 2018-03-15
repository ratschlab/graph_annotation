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
const char *shuf_seed = std::getenv("SHUF_SEED");
size_t seed = 0;

void serialize_vector_(std::ostream &out, const std::vector<cpp_int> &nums, size_t n_cols = 0) {
    n_cols += 64;
    for (auto it = nums.begin(); it != nums.end(); ++it) {
        size_t size = annotate::serialize(out, *it) + 8;
        if (size * 8 < n_cols) {
            char *zeros = (char*)malloc((n_cols + 7 - (size * 8)) >> 3);
            memset(zeros, 0, (n_cols + 7 - (size * 8)) >> 3);
            out.write(zeros, (n_cols + 7 - (size * 8)) >> 3);
            free(zeros);
        }
    }
}

uint64_t deserializeNumber(std::istream &in) {
    uint64_t n = 0;
    for (size_t i = 0; i < sizeof(n); ++i) {
        n = (n << 8) | in.get();
    }
    return n;
}

std::vector<cpp_int> read_comma_file(std::istream &in, size_t maxcount = -1llu) {
    std::vector<cpp_int> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    std::string line, digit;
    while (nums.size() != maxcount && std::getline(in, line)) {
        std::istringstream sin(line);
        nums.emplace_back(0);
        while (std::getline(sin, digit, ',')) {
            annotate::bit_set(nums.back(), std::stoi(digit));
            set_bits++;
        }
    }
    return nums;
}

std::vector<std::set<size_t>> read_comma_file_set(std::istream &in, size_t maxcount = -1llu) {
    std::vector<std::set<size_t>> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    std::string line, digit;
    while (nums.size() != maxcount && std::getline(in, line)) {
        std::istringstream sin(line);
        nums.emplace_back();
        while (std::getline(sin, digit, ',')) {
            annotate::bit_set(nums.back(), std::stoi(digit));
            set_bits++;
        }
    }
    return nums;
}

std::vector<std::vector<size_t>> read_comma_file_vector(std::istream &in, size_t maxcount = -1llu) {
    std::vector<std::vector<size_t>> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    std::string line, digit;
    while (nums.size() != maxcount && std::getline(in, line)) {
        std::istringstream sin(line);
        nums.emplace_back();
        while (std::getline(sin, digit, ',')) {
            annotate::bit_set(nums.back(), std::stoi(digit));
            set_bits++;
        }
    }
    return nums;
}


std::vector<cpp_int> read_strmap_file(
        std::istream &in,
        std::unordered_map<std::string, std::set<size_t>> &string_map,
        //std::unordered_map<std::string, std::set<size_t>>::iterator &strmap_it,
        size_t maxcount = -1llu,
        size_t memlim = (-1llu >> 30)) {
    //read from serialized index set
    std::vector<cpp_int> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    auto strmap_it = string_map.begin();
    while (nums.size() != maxcount && strmap_it != string_map.end()) {
        nums.emplace_back(0);
        for (auto &index : strmap_it->second) {
            annotate::bit_set(nums.back(), index);
            //set_bits++;
        }
        set_bits += strmap_it->second.size();
        strmap_it = string_map.erase(strmap_it);
        if (getCurrentRSS() > (memlim << 30))
            break;
        //++strmap_it;
    }
    return nums;
}

std::vector<std::set<size_t>> read_strmap_file_set(
        std::istream &in,
        std::unordered_map<std::string, std::set<size_t>> &string_map,
        size_t maxcount = -1llu) {
    std::vector<std::set<size_t>> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    auto strmap_it = string_map.begin();
    while (nums.size() != maxcount && strmap_it != string_map.end()) {
        nums.push_back(strmap_it->second);
        set_bits += strmap_it->second.size();
        strmap_it = string_map.erase(strmap_it);
    }
    return nums;
}

std::vector<std::vector<size_t>> read_strmap_file_vector(
        std::istream &in,
        std::unordered_map<std::string, std::set<size_t>> &string_map,
        size_t maxcount = -1llu) {
    std::vector<std::vector<size_t>> nums;
    if (maxcount < -1llu)
        nums.reserve(maxcount);
    auto strmap_it = string_map.begin();
    while (nums.size() != maxcount && strmap_it != string_map.end()) {
        std::vector<size_t> indices;
        indices.reserve(strmap_it->second.size());
        indices.insert(indices.end(), strmap_it->second.begin(), strmap_it->second.end());
        nums.push_back(indices);
        set_bits += strmap_it->second.size();
        strmap_it = string_map.erase(strmap_it);
    }
    return nums;
}



std::vector<cpp_int> read_raw_file(std::istream &in, size_t &num_rows, size_t maxcount = -1llu) {
    std::vector<cpp_int> nums;
    nums.reserve(std::min(num_rows, maxcount));
    while (nums.size() != maxcount && num_rows) {
        //nums.reserve(num_rows);
        size_t size = deserializeNumber(in);
        std::vector<uint64_t> row(size);
        for (auto it = row.begin(); it != row.end(); ++it) {
            *it = deserializeNumber(in);
            set_bits += __builtin_popcountll(*it);
        }
        if (shuf_seed) {
            //shuffle
            std::srand(seed);
            std::random_shuffle(
                    reinterpret_cast<uint8_t*>(row.data()),
                    reinterpret_cast<uint8_t*>(row.data() + row.size()),
                    [&](int i){ return std::rand()%i; });
        }
        total_bits += row.size() * 64;
        nums.emplace_back(0);
        mpz_import(nums.back().backend().data(), row.size(), -1, sizeof(row[0]), 0, 0, &row[0]);
        num_rows--;
    }
    return nums;
}

cpp_int to_num(const std::vector<size_t> &a) {
    cpp_int num = 0;
    for (auto &index : a) {
        annotate::bit_set(num, index);
    }
    return num;
}
/*
cpp_int to_num(const std::set<size_t> &a) {
    cpp_int num = 0;
    for (auto &index : a) {
        annotate::bit_set(num, index);
    }
    return num;
}
*/

int main(int argc, char** argv) {
    if (argc <= 2) {
        std::cerr << "ERROR: please pass in at least two files: an input and an output." << std::endl;
        exit(1);
    }
    std::string line, digit;
    std::vector<cpp_int> nums_ref;

    //environment variable arguments
    const char *step_char = std::getenv("STEP");
    size_t step = step_char ? atol(step_char) : -1llu;

    const char *test = std::getenv("TEST");

    const char *njobs = std::getenv("NJOBS");
    size_t n_jobs = 1;
    if (njobs) {
        n_jobs = atoi(njobs);
    }

    if (shuf_seed) {
        seed = atoi(shuf_seed);
    }
    //const char *dump_raw = std::getenv("DUMP");
    /*
    std::ofstream dout;
    if (dump_raw) {
        dump_cols = atoi(dump_raw);
        dout.open(std::string(argv[argc - 1]) + ".raw");
    }
    *///size_t dump_cols = 0;

    const char *read_comma = std::getenv("COMMA");

    const char *strmap = std::getenv("MAP");

    const char *memlim = std::getenv("MEM");

    const char *indexset = std::getenv("INDEXSET");
    size_t mem_lim = -1llu >> 30;
    if (memlim) {
        mem_lim = atoi(memlim);
    }

    const char *maxtest = std::getenv("MAXTEST");
    size_t max_test = 0;
    if (maxtest) {
        max_test = atol(maxtest);
    }

    //std::vector<annotate::WaveletTrie*> wtrs;
    std::vector<std::future<annotate::WaveletTrie*>> wtrs;

    double runtime = 0;
    double readtime = 0;
    //double dumptime = 0;
    size_t num_rows = 0;
    for (int f = 1; f < argc - 1; ++f) {
        std::ifstream fin(argv[f]);
        if (!fin.good()) {
            std::cerr << "WARNING: file " << argv[f] << " bad." << std::endl;
            exit(1);
        }
        if (!read_comma && !strmap) {
            num_rows = deserializeNumber(fin);
            if (n_jobs > 1) {
                step = std::min(num_rows / n_jobs + 1, step);
            }
        }
        std::unordered_map<std::string, std::set<size_t>> string_map;
        if (strmap) {
            boost::archive::binary_iarchive iarch(fin);
            size_t num_elements;
            std::cout << "Loading input" << std::endl;
            iarch & string_map;
            iarch & num_elements;
            if (n_jobs > 1) {
                step = std::min(string_map.size() / n_jobs + 1, step);
            }
        }
        //auto strmap_it = string_map.begin();

        std::cout << "Compressing:\t" << argv[f] << std::endl;
        std::cout << "Step size:\t" << step << std::endl;
#ifndef NPRINT
        auto policy = std::launch::deferred;
#else
        auto policy = std::launch::async;
#endif
        if (indexset) {
            while (true) {
                std::vector<std::vector<size_t>> *nums;
                if (strmap) {
                    nums = new std::vector<std::vector<size_t>>(read_strmap_file_vector(fin, string_map, step));
                } else if (read_comma) {
                    nums = new std::vector<std::vector<size_t>>(read_comma_file_vector(fin, step));
                } else {
                    std::cerr << "Only precise annotator supported\n";
                    exit(1);
                }
                if (!nums->size()) {
                    delete nums;
                    break;
                }
                if (test != NULL && nums_ref.size() < max_test) {
                    nums_ref.reserve(nums_ref.size() + nums->size());
                    std::transform(nums->begin(), nums->end(), std::back_inserter(nums_ref), to_num);
                }
                wtrs.push_back(std::async(policy, [=]() {
                    std::cout << "s" << std::flush;
                    auto wtr = new annotate::WaveletTrie(nums->begin(), nums->end());
                    delete nums;
                    std::cout << "." << std::flush;
                    return wtr;
                }));
            }
        } else {
        while (true) {
            std::vector<cpp_int> *nums;
            if (read_comma) {
                //read from text index list
                nums = new std::vector<cpp_int>(read_comma_file(fin, step));
            } else if (strmap) {
                //read from serialized index set
                nums = new std::vector<cpp_int>(read_strmap_file(fin, string_map, step, mem_lim));
                //nums = read_strmap_file(fin, string_map, strmap_it, step);
            } else {
                nums = new std::vector<cpp_int>(read_raw_file(fin, num_rows, step));
            }
            if (!nums->size()) {
                delete nums;
                break;
            }
            if (test != NULL) {
                nums_ref.reserve(nums_ref.size() + nums->size());
                nums_ref.insert(nums_ref.end(), nums->begin(), nums->end());
            }
            //wtrs.push_back(new annotate::WaveletTrie(nums));
            wtrs.push_back(std::async(policy, [=]() {
                std::cout << "s" << std::flush;
                auto wtr = new annotate::WaveletTrie(nums->begin(), nums->end());
                delete nums;
                std::cout << "." << std::flush;
                return wtr;
            }));
            if (wtrs.size() % n_jobs == 0 || getCurrentRSS() > (mem_lim << 30)) {
                //std::cout << "j" << std::flush;
                for (auto it = wtrs.rbegin(); it != wtrs.rbegin() + n_jobs; ++it) {
                    it->wait();
                }
                if (getCurrentRSS() > (mem_lim << 30)) {
                    std::cerr << "ERROR: memory limit set too low" << std::endl;
                    exit(1);
                }
            }
            //std::cout << "." << std::flush;
        }
        }
        //std::cout << std::endl;
        fin.close();
    }
    std::cout << std::endl;

    //if (wtrs.size() == 0)
    //    return 0;

    //if (n_jobs > 1)
        //omp_set_num_threads(n_jobs);

    //auto *wtr = wtrs[0];
    annotate::WaveletTrie *wtr = new annotate::WaveletTrie();
    if (wtrs.size()) {
        //wtr = wtrs[0].get();
        std::cout << "Merging" << std::endl;
        for (auto it = wtrs.begin(); it != wtrs.end(); ++it) {
            auto *curwtr = it->get();
            wtr->insert(*curwtr);
            delete curwtr;
            //wtr->insert(**it);
            //wtrs[0]->insert(**it);
            //delete *it;
            std::cout << "." << std::flush;
        }
    }
    std::cout << std::endl;

    std::cout << "Times:" << std::endl;
    std::cout << "Reading:\t" << readtime << std::endl;
    std::cout << "Compressing:\t" << runtime << std::endl;

    if (test != NULL) {
        std::cout << "Decompressing:\t" << std::flush;
        assert(wtr->size() == nums_ref.size());
#ifndef NPRINT
        wtr->print();
#endif
        Timer decomp_timer;
        decomp_timer.reset();
        for (size_t i = 0; i < nums_ref.size(); ++i) {
            if (nums_ref.at(i) != wtr->at(i)) {
                std::cerr << "Fail at " << i << "\n";
                std::cerr << nums_ref.at(i) << "\n" << wtr->at(i) << "\n";
                exit(1);
            }
        }
        double decom_time = decomp_timer.elapsed();
        std::cout << "Check time:\t" << decom_time << "\n";
        std::cout << "Time per edge:\t" << decom_time / nums_ref.size() << "\n";
        std::cout << std::endl;
    }
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
        if (shuf_seed) {
            std::cout << "Shuffle seed:\t" << seed << std::endl;
        }
        fout.close();
        if (test != NULL) {
            std::cout << "Testing deserialization:\t" << std::flush;
            std::ifstream fin(argv[argc - 1]);
            annotate::WaveletTrie wtr_load;
            wtr_load.load(fin);
            fin.close();
            if (*wtr != wtr_load) {
                std::cerr << "Serialization/deserialization failed" << std::endl;
                exit(1);
            }
            std::cout << std::endl;
        }
        delete wtr;
    }
    /*
    if (dump_raw) {
        dout.close();
        std::cout << "Raw dump:\t" << dumptime << std::endl;
    }
    */
    std::cout << "Done\n";
    return 0;
}


