#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

namespace serialization {
    uint64_t serializeNumber(std::ostream &out, uint64_t const n);

    uint64_t serializeNumberVector(std::ostream &out, std::vector<uint64_t> const &v);

    uint64_t loadNumber(std::istream &in);

    std::vector<uint64_t> loadNumberVector(std::istream &in);

    uint64_t serializeString(std::ostream &out, const std::string &s);

    std::string loadString(std::istream &in);

    uint64_t serializeStringMap(
            std::ostream &out,
            const std::unordered_map<std::string, size_t> &umap);

    std::unordered_map<std::string, size_t> loadStringMap(std::istream &in);

} // namespace serialization

#endif // __SERIALIZATION_HPP__
