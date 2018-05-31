#ifndef __WAVELET_TRIE_ANNOTATOR__
#define __WAVELET_TRIE_ANNOTATOR__

#include <fstream>

#include "wavelet_trie.hpp"
#include "dbg_bloom_annotator.hpp"
#include "serialization.hpp"


namespace annotate {

class WaveletTrieAnnotator {
  public:
    WaveletTrieAnnotator(const hash_annotate::DeBruijnGraphWrapper &graph, size_t p = 1)
        : graph_(graph),
          wt_(new WaveletTrie(p)),
          num_columns_(0) {}
    
    WaveletTrieAnnotator(const hash_annotate::PreciseHashAnnotator &precise,
                         const hash_annotate::DeBruijnGraphWrapper &graph,
                         size_t p = 1)
        : graph_(graph),
          num_columns_(precise.num_columns()) {
        auto annots = extract_raw_annots(precise);
        wt_ = new WaveletTrie(annots, p);
    }

    ~WaveletTrieAnnotator() { delete wt_; }

    std::vector<uint64_t> annotate_edge(hash_annotate::DeBruijnGraphWrapper::edge_index i) const {
        auto vect = wt_->at(i);
        size_t a = mpz_size(vect.backend().data());
        std::vector<uint64_t> ret_vect((a + 63) >> 3);
        mpz_export(ret_vect.data(), &a, -1, sizeof(uint64_t), 0, 0, vect.backend().data());
#ifndef NDEBUG
        if (((num_columns_ + 63) >> 6) < ret_vect.size()) {
            for (auto it = ret_vect.begin() + a; it != ret_vect.end(); ++it) {
                assert(!*it);
            }
        }
#endif
        ret_vect.resize((num_columns_ + 63) >> 6);
        //memcpy(reinterpret_cast<char*>(ret_vect.data()), l_int_raw, a);
        //free(l_int_raw);
        return ret_vect;
    }

    uint64_t serialize(std::ostream &out) const {
        return wt_->serialize(out)
             + serialization::serializeNumber(out, num_columns_);
    }
    void load(std::istream &in) {
        wt_->load(in);
        num_columns_ = serialization::loadNumber(in);
    }

    uint64_t serialize(const std::string &filename) const {
        std::ofstream out(filename);
        return serialize(out);
    }

    void load(const std::string &filename) {
        std::ifstream in(filename);
        load(in);
    }

    size_t size() const { return wt_->size(); }

  private:
    const hash_annotate::DeBruijnGraphWrapper &graph_;
    WaveletTrie *wt_;
    size_t num_columns_;

    std::vector<cpp_int> extract_raw_annots(const hash_annotate::PreciseHashAnnotator &precise) {
        std::vector<cpp_int> annots;
        auto permut_map = precise.compute_permutation_map();
        // order determined by hash map
        for (auto &kmer : precise) {
            annots.emplace_back(0);
            for (auto &index : kmer.second) {
                bit_set(annots.back(), permut_map.size() ? permut_map[index] : index);
            }
        }

        // alt: order determined by graph
        /*
        for (size_t i = 0; i < precise.size(); ++i) {
            annots.emplace_back(0);
            auto vect = precise.permute_indices(precise.annotate_edge(i), permut_map);
            mpz_import(annots.back().backend().data(), vect.size(), -1, sizeof(vect[0]), 0, 0, &vect[0]);
        }
        */
        return annots;
    }
};
}; // namespace annotate

#endif // __WAVELET_TRIE_ANNOTATOR__
