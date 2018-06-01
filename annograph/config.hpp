#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>


class Config {
  public:
    Config(int argc, const char *argv[]);

    bool verbose = false;
    bool reverse = false;
    bool fasta_anno = false;
    bool wavelet_trie = false;

    unsigned int k = 3;
    unsigned int distance = 0;
    unsigned int alignment_length = 0;
    unsigned int bloom_num_hash_functions = 0;
    unsigned int bloom_test_num_kmers = 0;
    unsigned int p = 1;

    double bloom_fpp = -1;
    double bloom_bits_per_edge = -1;
    double discovery_fraction = 1.0;

    std::vector<std::string> fname;
    std::string outfbase;
    std::string infbase;
    std::string dbpath;
    std::string refpath;
    std::string fasta_header_delimiter;

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD,
        MAP,
    };
    IdentityType identity = NO_IDENTITY;

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
