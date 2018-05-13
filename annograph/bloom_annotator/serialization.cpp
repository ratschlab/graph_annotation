#include "serialization.hpp"

namespace serialization {
    uint64_t serializeNumber(std::ostream &out, uint64_t const n) {
        out.put(static_cast<char>((n >> (7 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (6 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (5 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (4 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (3 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (2 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (1 * 8)) & 0xFF));
        out.put(static_cast<char>((n >> (0 * 8)) & 0xFF));

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

    uint64_t loadNumber(std::istream &in) {
        uint64_t number = 0;
        number |= static_cast<uint64_t>(in.get()) << (7 * 8);
        number |= static_cast<uint64_t>(in.get()) << (6 * 8);
        number |= static_cast<uint64_t>(in.get()) << (5 * 8);
        number |= static_cast<uint64_t>(in.get()) << (4 * 8);
        number |= static_cast<uint64_t>(in.get()) << (3 * 8);
        number |= static_cast<uint64_t>(in.get()) << (2 * 8);
        number |= static_cast<uint64_t>(in.get()) << (1 * 8);
        number |= static_cast<uint64_t>(in.get()) << (0 * 8);
        return number;
    }

    std::vector<uint64_t> loadNumberVector(std::istream &in) {
        uint64_t s = loadNumber(in);
        std::vector<uint64_t> vect(s);
        for (auto &n : vect) {
            n = loadNumber(in);
        }
        return vect;
    }

    uint64_t serializeString(std::ostream &out, const std::string &s) {
        out.write(&s[0], s.length() + 1);
        return s.length() + 1;
    }

    std::string loadString(std::istream &in) {
        std::vector<char> s;
        do { s.push_back(in.get()); } while (s.back());
        return std::string(s.begin(), s.end() - 1);
    }

    uint64_t serializeStringMap(
            std::ostream &out,
            const std::unordered_map<std::string, size_t> &umap) {
        size_t total_size = 0;
        total_size += serializeNumber(out, umap.size());
        for (const auto &item : umap) {
            total_size += serializeString(out, item.first);
            total_size += serializeNumber(out, item.second);
        }
        return total_size;
    }

    std::unordered_map<std::string, size_t> loadStringMap(std::istream &in) {
        std::unordered_map<std::string, size_t> umap;
        size_t total_size = loadNumber(in);
        while (total_size--) {
            std::string key = loadString(in);
            size_t value = loadNumber(in);
            umap[key] = value;
        }
        return umap;
    }

} // namespace serialization

