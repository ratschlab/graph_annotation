#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <unordered_set>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/impl/basic_binary_oprimitive.ipp>
#include <boost/archive/impl/basic_binary_iprimitive.ipp>

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
    if (argc < 2) {
        std::cerr << argc << std::endl;
        std::cerr << "Too few arguments" << std::endl;
        exit(1);
    }
    std::ifstream fin(argv[1]);
    if (std::getenv("RAW")) {
        size_t num_rows = deserializeNumber(fin);
        std::cout << "Num elements:\t" << num_rows << std::endl;
        std::unordered_set<std::string> elements;
        while (num_rows) {
            size_t size = deserializeNumber(fin);
            std::string curstring;
            size_t index;
            while (size) {
                index = deserializeNumber(fin);
                curstring += std::to_string(index) + ",";
                size--;
            }
            elements.insert(curstring);
            num_rows--;
        }
        std::cout << "Num uniq elements:\t" << elements.size() << std::endl;
        fin.close();
        return 0;
    }
    boost::archive::binary_iarchive iarch(fin);
    std::unordered_map<std::string, std::set<size_t>> string_map;
    size_t num_elements;
    std::cout << "Loading input" << std::endl;
    iarch & string_map;
    iarch & num_elements;
    fin.close();
    std::cout << "Num elements:\t" << string_map.size() << std::endl;
    std::unordered_set<std::string> elements;
    std::unordered_set<std::string> kmers;
    for (auto &map : string_map) {
        kmers.insert(map.first.substr(0,map.first.length()-1));
        std::string curstring;
        for (auto &index : map.second) {
            curstring += std::to_string(index) + ",";
        }
        elements.insert(curstring);
    }
    std::cout << "Num nodes:\t" << kmers.size() << std::endl;
    std::cout << "Num uniq elements:\t" << elements.size() << std::endl;
    return 0;
}
