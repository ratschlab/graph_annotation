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
          wt_(p),
          num_columns_(0) {}
    
    WaveletTrieAnnotator(const hash_annotate::PreciseHashAnnotator &precise,
                         const hash_annotate::DeBruijnGraphWrapper &graph,
                         size_t p = 1)
        : graph_(graph),
          num_columns_(precise.num_columns()) {
        if (p == 1) {
            wt_ = annotate::WaveletTrie(extract_index_set(precise), p);
            //wt_ = annotate::WaveletTrie(extract_raw_annots(precise), p);
        } else {
            utils::ThreadPool thread_queue(p);
            size_t step = (precise.size() + p - 1) / p;
            std::vector<std::future<annotate::WaveletTrie>> wtrs;
            std::vector<std::vector<size_t>> indices;
            indices.reserve(step);

            for (size_t i = 0; i < precise.size(); ++i) {
                auto kmer_indices = precise.annotate_edge_indices(i);
                indices.emplace_back(kmer_indices.begin(), kmer_indices.end());
                if (indices.size() == step) {
                    wtrs.emplace_back(
                        thread_queue.enqueue([](std::vector<std::vector<size_t>> indices) {
                            return annotate::WaveletTrie(std::move(indices));
                        }, indices)
                    );
                    indices.clear();
                }
            }
            if (indices.size()) {
                wtrs.emplace_back(
                    thread_queue.enqueue([](std::vector<std::vector<size_t>> &indices) {
                        return annotate::WaveletTrie(std::move(indices));
                    }, indices)
                );
            }

            for (auto &wtr : wtrs) {
                wt_.insert(wtr.get());
            }
        }
    }

    //TODO rewrite as load_from_precise
    WaveletTrieAnnotator(std::istream &in,
                         const hash_annotate::DeBruijnGraphWrapper &graph,
                         size_t p = 1)
        : graph_(graph),
          wt_(p) {
        hash_annotate::PreciseHashAnnotator precise(graph_);
        precise.load(in);
        num_columns_ = precise.num_columns();
        wt_ = annotate::WaveletTrie(extract_index_set(precise), p);
        /*
        if (false && p == 1) {
            hash_annotate::PreciseHashAnnotator precise(graph);
            precise.load(in);
            num_columns_ = precise.num_columns();
            wt_ = annotate::WaveletTrie(extract_index_set(precise), p);
        } else {
            // TODO: this depends too much PreciseHashAnnotator serialization
            //       is implemented.
            size_t prefix_size = serialization::loadNumber(in);
            std::set<uint64_t> prefix_indices;
            while (prefix_size--) {
                prefix_indices.insert(serialization::loadNumber(in));
            }

            num_columns_ = serialization::loadNumber(in);

            auto permut_map = hash_annotate::PreciseHashAnnotator::compute_permutation_map(
                    num_columns_,
                    prefix_indices);

            size_t graph_size = serialization::loadNumber(in);

            utils::ThreadPool thread_queue(p);
            size_t step = (graph_size + p - 1) / p;
            std::vector<std::future<annotate::WaveletTrie>> wtrs;
            std::vector<std::vector<size_t>> indices;
            indices.reserve(step);

            while (graph_size--) {
                indices.emplace_back();
                std::set<size_t> nums;
                size_t num_size = serialization::loadNumber(in);
                while (num_size--) {
                    indices.back().emplace_back(prefix_indices.size()
                            ? permut_map[serialization::loadNumber(in)]
                            : serialization::loadNumber(in));
                }

                //TODO: optional?
                std::sort(indices.back().begin(), indices.back().end());

                auto dummy = serialization::loadString(in);
                if (indices.size() == step) {
                    wtrs.emplace_back(
                        thread_queue.enqueue([](std::vector<std::vector<size_t>> indices) {
                            return annotate::WaveletTrie(std::move(indices));
                        }, indices)
                    );
                    indices.clear();
                }
            }
            if (indices.size()) {
                wtrs.emplace_back(
                    thread_queue.enqueue([](std::vector<std::vector<size_t>> &indices) {
                        return annotate::WaveletTrie(std::move(indices));
                    }, indices)
                );
            }

            for (auto &wtr : wtrs) {
                wt_.insert(wtr.get());
            }
        }
        */
    }

    ~WaveletTrieAnnotator() {}

    std::vector<uint64_t> annotate_edge(hash_annotate::DeBruijnGraphWrapper::edge_index i) const {
        auto vect = wt_.at(i);
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
        //return serialization::serializeNumber(out, num_columns_)
        //     + wt_.serialize(out);
        return wt_.serialize(out) + serialization::serializeNumber(out, num_columns_);
    }
    uint64_t serialize(const std::string &filename) const {
        std::ofstream out(filename);
        return serialize(out);
    }

    void load(std::istream &in) {
        //num_columns_ = serialization::loadNumber(in);
        //wt_.load(in);
        wt_.load(in);
        num_columns_ = serialization::loadNumber(in);
    }
    void load(const std::string &filename) {
        std::ifstream in(filename);
        load(in);
    }

    size_t size() const { return wt_.size(); }

    size_t num_columns() const { return num_columns_; }

    bool operator==(const WaveletTrieAnnotator &that) const {
        return num_columns_ == that.num_columns_
            && wt_ == that.wt_;
    }

  private:
    const hash_annotate::DeBruijnGraphWrapper &graph_;
    WaveletTrie wt_;
    size_t num_columns_;

    std::vector<cpp_int> extract_raw_annots(const hash_annotate::PreciseHashAnnotator &precise) {
        std::vector<cpp_int> annots;
        annots.reserve(precise.size());
        auto permut_map = precise.compute_permutation_map();

        /*
        // order determined by hash map
        for (auto &kmer : precise) {
            for (auto &index : kmer.second) {
                bit_set(annots.back(), permut_map.size() ? permut_map[index] : index);
            }
        }
        */

        // alt: order determined by graph
        for (size_t i = 0; i < precise.size(); ++i) {
            annots.emplace_back(0);
            auto vect = precise.permute_indices(precise.annotate_edge(i), permut_map);
            mpz_import(annots.back().backend().data(), vect.size(), -1, sizeof(vect[0]), 0, 0, &vect[0]);
        }
        return annots;
    }

    std::vector<std::vector<size_t>> extract_index_set(const hash_annotate::PreciseHashAnnotator &precise) {
        std::vector<std::vector<size_t>> indices;
        indices.reserve(precise.size());
        auto permut_map = precise.compute_permutation_map();
        //for (auto &kmer : precise) {
        for (size_t j = 0; j < precise.size(); ++j) {
            auto kmer_indices = precise.annotate_edge_indices(j);
            if (permut_map.empty()) {
                indices.emplace_back(kmer_indices.begin(), kmer_indices.end());
            } else {
                indices.emplace_back();
                for (auto i : kmer_indices) {
                    indices.back().emplace_back(permut_map[i]);
                }
                std::sort(indices.back().begin(), indices.back().end());
            }
        }
        return indices;
    }

};
}; // namespace annotate

#endif // __WAVELET_TRIE_ANNOTATOR__
