#include "config.hpp"

#include <cstring>
#include <iostream>


Config::Config(int argc, const char *argv[]) {
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(argv[0]);
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "build")) {
        identity = BUILD;
    } else if (!strcmp(argv[1], "update")) {
        identity = UPDATE;
    } else if (!strcmp(argv[1], "map")) {
        identity = MAP;
    } else if (!strcmp(argv[1], "permutation")) {
        identity = PERMUTATION;
    } else if (!strcmp(argv[1], "stats")) {
        identity = STATS;
    } else if (!strcmp(argv[1], "query")) {
        identity = QUERY;
    }
    // provide help screen for chosen identity
    if (argc == 2) {
        print_usage(argv[0], identity);
        exit(-1);
    }

    // parse remaining command line items
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reverse")) {
            reverse = true;
        } else if (!strcmp(argv[i], "--fasta-anno")) {
            fasta_anno = true;
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            p = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--wavelet-trie")) {
            wavelet_trie = true;
        } else if (!strcmp(argv[i], "--bloom-false-pos-prob")) {
            bloom_fpp = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--bloom-bits-per-edge")) {
            bloom_bits_per_edge = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--discovery-fraction")) {
            discovery_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--bloom-hash-functions")) {
            bloom_num_hash_functions = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--bloom-test-num-kmers")) {
            bloom_test_num_kmers = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--num-permutations")) {
            num_permutations = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--reference")) {
            refpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--fasta-header-delimiter")) {
            fasta_header_delimiter = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0], identity);
            exit(0);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\nERROR: Unknown option %s\n\n", argv[i]);
            print_usage(argv[0], identity);
            exit(-1);
        } else {
            fname.push_back(argv[i]);
        }
    }

    /*
    if (!fname.size()) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.size())
                fname.push_back(line);
        }
    }
    */

    bool print_usage_and_exit = false;

    if (identity == PERMUTATION && infbase.empty())
        print_usage_and_exit = true;

    if (identity == MAP && !fname.size())
        print_usage_and_exit = true;

    if (identity == UPDATE && !infbase.size())
        print_usage_and_exit = true;

    if (!fname.size() && infbase.empty())
        print_usage_and_exit = true;

    if (fasta_header_delimiter.size() > 1) {
        std::cerr << "FASTA header delimiter must be at most one character" << std::endl;
        print_usage_and_exit = true;
    }

    // if misused, provide help screen for chosen identity and exit
    if (print_usage_and_exit) {
        print_usage(argv[0], identity);
        exit(-1);
    }
}

void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    fprintf(stderr, "Metagenome graph representation -- Version 0.1\n\n");

    switch (identity) {
        case NO_IDENTITY: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available commands:\n\n");

            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tupdate\t\textend graph from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tmap\t\tannotate a k-mer\n\n");
            return;
        }
        case BUILD: {
            fprintf(stderr, "Usage: %s build [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t   --reference [STR] \t\t\tbasename of reference sequence []\n");
            fprintf(stderr, "\t   --fasta-header-delimiter [STR] \theader delimiter (for setting multiple annotations) []\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\t\t\tbasename of output file []\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \t\t\tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\t\tnumber of threads to use for wavelet trie compression [1]\n");
            fprintf(stderr, "\t   --wavelet-trie \t\t\tconstruct wavelet trie [off]\n");
            fprintf(stderr, "\t   --bloom-false-pos-prob [FLOAT] \tFalse positive probability in bloom filter [-1]\n");
            fprintf(stderr, "\t   --bloom-bits-per-edge [FLOAT] \tBits per edge used in bloom filter annotator [0.4]\n");
            fprintf(stderr, "\t   --bloom-hash-functions [INT] \tNumber of hash functions used in bloom filter [off]\n");
            fprintf(stderr, "\t   --bloom-test-num-kmers \t\tEstimate false positive rate for every n k-mers [0]\n");
            fprintf(stderr, "\t-r --reverse \t\t\t\tadd reverse complement reads [off]\n");
        } break;
        case UPDATE: {
            fprintf(stderr, "Usage: %s update [options] -i <graph_basename> FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for update:\n");
            fprintf(stderr, "\t   --reference [STR] \t\t\tbasename of reference sequence []\n");
            fprintf(stderr, "\t   --fasta-header-delimiter [STR] \theader delimiter (for setting multiple annotations) []\n");
            fprintf(stderr, "\t-i --infile-base [STR] \t\t\tinput colored graph basename\n");
            fprintf(stderr, "\t   --wavelet-trie \t\t\tuse wavelet trie for annotation [off]\n"
                            "\t                  \t\t\t                         default: Bloom filter\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\t\t\tbasename of output file []\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\t\tnumber of threads to use for wavelet trie compression [1]\n");
            fprintf(stderr, "\t-r --reverse \t\t\t\tadd reverse complement reads [off]\n");
        } break;
        case MAP: {
            fprintf(stderr, "Usage: %s map [options] -i <graph_basename> k-mer [[k-mer] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for map:\n");
            fprintf(stderr, "\t   --wavelet-trie \tuse wavelet trie for annotation [off]\n"
                            "\t                  \t                         default: Bloom filter\n");
            fprintf(stderr, "\t-i --infile-base [STR] \tinput colored graph basename\n");
            // fprintf(stderr, "\t-p --parallel [INT] \tnumber of threads to use for wavelet trie compression [1]\n");
        } break;
        case PERMUTATION: {
            fprintf(stderr, "Usage: %s permutation [options] -i <graph_basename>\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for permutation:\n");
            fprintf(stderr, "\t   --num-permutations \t\tnumber of column index permutations to test [0]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tnumber of threads (one permutation per thread) [1]\n");
            fprintf(stderr, "\t-i --infile-base [STR] \tinput colored graph basename\n");
        } break;
        case QUERY: {
        } break;
        case STATS: {
        } break;
    }

    fprintf(stderr, "\n\tGeneral options:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}
