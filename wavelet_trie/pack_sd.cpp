#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include <sdsl/wavelet_trees.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>
#define BVSIZE 1024

using namespace std;

uint64_t deserializeNumber(std::istream &in) {
    uint64_t n = 0;
    for (size_t i = 0; i < sizeof(n); ++i) {
        if (in.eof()) {
            std::cerr << "ERROR: reached end of file" << std::endl;
            exit(1);
        }
        n = (n << 8) | in.get();
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << argc << std::endl;
        std::cerr << "Too few arguments" << std::endl;
        exit(1);
    }
    std::ifstream fin(argv[1]);
    std::ofstream fout(argv[2]);
    std::ofstream fcout(argv[3]);
    if (std::getenv("MAP")) {
        boost::archive::binary_iarchive iarch(fin);
        boost::archive::binary_oarchive oarch(fcout);
        std::unordered_map<std::string, std::set<size_t>> string_map;
        size_t num_elements;
        std::cout << "Loading input" << std::endl;
        iarch & string_map;
        iarch & num_elements;
        fin.close();
        size_t set_bits = 0;
        for (auto &curstring : string_map) {
            set_bits += curstring.second.size();
        }
        std::cout << "Rows:\t" << string_map.size() << std::endl;
        std::cout << "Num cols:\t" << num_elements << std::endl;
        std::cout << "Total bits:\t" << string_map.size() * num_elements << std::endl;
        std::cout << "Set bits:\t" << set_bits << std::endl;
        sdsl::sd_vector_builder builder(string_map.size() * num_elements, set_bits);
        size_t counter = 0;
        for (auto &curstring : string_map) {
            std::vector<uint64_t> raw_annot((num_elements + 63) >> 6);
            for (auto &index : curstring.second) {
                builder.set(counter + index);
                raw_annot[index >> 6] |= 1llu << (index % 64);
            }
            counter += num_elements;
            oarch & raw_annot;
        }
        sdsl::sd_vector<>(builder).serialize(fout);
        fout.close();
        fcout.close();
        return 0;
    }
    size_t num_rows = deserializeNumber(fin);
    std::cout << "Rows:\t" << num_rows << std::endl;
    size_t set_bits;
    const char* setbits = std::getenv("SETBITS");
    if (setbits) {
        set_bits = atol(setbits);
        std::cout << "Row size:\t" << BVSIZE << std::endl;
        std::cout << "Set bits:\t" << set_bits << std::endl;
        sdsl::sd_vector_builder builder(num_rows * BVSIZE, set_bits);
        size_t counter = 0;
        while (num_rows--) {
            size_t size = deserializeNumber(fin);
            if (size > (BVSIZE >> 6)) {
                std::cerr << "More that " << BVSIZE << "bits" << std::endl;
            }
            for (size_t i = 0; i < (BVSIZE >> 6); ++i) {
                size_t limb;
                if (i < size) {
                    limb = deserializeNumber(fin);
                    size_t j = 0;
                    for (size_t jj = limb; jj; jj >>= 1) {
                        if (jj & 1) {
                            builder.set(counter + j);
                        }
                        j++;
                    }
                    counter += BVSIZE;
                } else {
                    limb = 0;
                }
                fcout.write(reinterpret_cast<const char*>(&limb), sizeof(limb));
            }
        }
        sdsl::sd_vector<>(builder).serialize(fout);
        fout.close();
        fcout.close();
        return 0;
    }
    sdsl::bit_vector bits(num_rows * BVSIZE);
    size_t counter = 0;
    while (num_rows--) {
        size_t size = deserializeNumber(fin);
        if (size > (BVSIZE >> 6)) {
            std::cerr << "More that " << BVSIZE << "bits" << std::endl;
        }
        for (size_t i = 0; i < (BVSIZE >> 6); ++i) {
            size_t limb;
            if (i < size) {
                limb = deserializeNumber(fin);
                bits.set_int(counter, limb);
                counter += BVSIZE;
            } else {
                limb = 0;
            }
            fcout.write(reinterpret_cast<const char*>(&limb), sizeof(limb));
        }
    }
    sdsl::sd_vector<>(bits).serialize(fout);
    fout.close();
    fcout.close();
    return 0;
}
