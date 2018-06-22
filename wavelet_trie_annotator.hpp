#ifndef __WAVELET_TRIE_ANNOTATOR__
#define __WAVELET_TRIE_ANNOTATOR__

#include <fstream>

#include "wavelet_trie.hpp"
#include "dbg_bloom_annotator.hpp"
#include "serialization.hpp"

namespace annotate {

class WaveletTrieAnnotator : public hash_annotate::PreciseAnnotator {
  public:
    WaveletTrieAnnotator(const hash_annotate::DeBruijnGraphWrapper &graph, size_t p = 1);
    
    WaveletTrieAnnotator(const hash_annotate::PreciseHashAnnotator &precise,
                         const hash_annotate::DeBruijnGraphWrapper &graph,
                         size_t p = 1,
                         std::unordered_map<pos_t, pos_t>&& permut_map = {});

    void load_from_precise_file(std::istream &in, size_t p = 1);

    ~WaveletTrieAnnotator() {}

    void add_sequence(const std::string &sequence,
                      hash_annotate::pos_t column = static_cast<hash_annotate::pos_t>(-1),
                      bool rooted = false);

    void add_column(const std::string &sequence, bool rooted = false);

    std::vector<uint64_t> annotate_edge(hash_annotate::DeBruijnGraphWrapper::edge_index i, bool permute = false) const;

    std::vector<uint64_t> annotation_from_kmer(const std::string &kmer, bool permute = false) const;

    uint64_t serialize(std::ostream &out) const;
    uint64_t serialize(const std::string &filename) const;

    bool load(std::istream &in);
    bool load(const std::string &filename);

    size_t size() const { return wt_.size(); }

    size_t num_columns() const { return num_columns_; }

    bool operator==(const WaveletTrieAnnotator &that) const {
        return num_columns_ == that.num_columns_
            && wt_ == that.wt_;
    }
    bool operator!=(const WaveletTrieAnnotator &that) const { return !(operator==(that)); }

    void print() const { wt_.print(); }

    std::tuple<size_t, size_t, size_t, size_t> stats() const {
        auto num_uniq_set_bits = wt_.stats();
        return std::make_tuple(size(), num_columns(), num_uniq_set_bits.second, num_uniq_set_bits.first);
    }

  private:
    const hash_annotate::DeBruijnGraphWrapper &graph_;
    WaveletTrie wt_;
    size_t num_columns_;
    std::unordered_map<pos_t, pos_t> permut_map_;

    std::vector<cpp_int> extract_raw_annots(const hash_annotate::PreciseHashAnnotator &precise);

    std::vector<std::vector<pos_t>> extract_index_set(
            const hash_annotate::PreciseHashAnnotator &precise,
            std::unordered_map<pos_t, pos_t>&& permut_map = {});

};
}; // namespace annotate

#endif // __WAVELET_TRIE_ANNOTATOR__
