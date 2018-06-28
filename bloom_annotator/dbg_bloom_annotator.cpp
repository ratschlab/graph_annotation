#include "dbg_bloom_annotator.hpp"

#include <fstream>
#include <cmath>
#include <map>
#include <unordered_map>

#include "../serialization.hpp"


namespace hash_annotate {

std::unordered_map<size_t, size_t> PreciseHashAnnotator::compute_permutation_map() const {
    std::unordered_map<size_t, size_t> index_map;
    size_t index_size = 0;
    for (auto i : prefix_indices_) {
        // indices are wrong when setting index_map[i] = index_map.size() instead
        index_map[i] = index_size++;
    }
    if (index_map.empty())
        return index_map;
    for (size_t i = 0; i < annotation_exact.size(); ++i) {
        index_size += index_map.emplace(i, index_size).second;
    }
    assert(index_size == num_columns());
    /*
    std::vector<pos_t> indices(annotation_exact.size());
    std::iota(indices.begin(), indices.end(), 0);
    auto it = indices.begin();
    for (auto i : prefix_indices_) {
        std::rotate(it, it + 1, indices.begin() + i + 1);
        ++it;
    }
    */
    return index_map;
}

std::unordered_map<size_t, size_t>
PreciseHashAnnotator::compute_permutation_map(size_t num_columns,
                                              const std::set<pos_t> &prefix_indices) {
    std::unordered_map<size_t, size_t> index_map;
    size_t index_size = 0;
    for (auto i : prefix_indices) {
        index_map[i] = index_size++;
    }
    if (index_map.empty())
        return index_map;
    for (size_t i = 0; i < num_columns; ++i) {
        index_size += index_map.emplace(i, index_size).second;
    }
    assert(index_size == num_columns);
    return index_map;
}

std::vector<uint64_t>
PreciseHashAnnotator::permute_indices(const std::vector<uint64_t> &a,
                                      const std::unordered_map<size_t, size_t> &index_map) const {
    if (index_map.empty())
        return a;
    std::vector<uint64_t> b(a.size());
    for (size_t i = 0; i < annotation_exact.size(); ++i) {
        if (test_bit(a, i))
            set_bit(b, index_map.find(i)->second);
    }
    return b;
}

uint64_t PreciseHashAnnotator::serialize(std::ostream &out) const {
    uint64_t written_bytes = 0;

    written_bytes += serialization::serializeNumber(out, prefix_indices_.size());
    for (auto i : prefix_indices_)
        written_bytes += serialization::serializeNumber(out, i);

    written_bytes += annotation_exact.serialize(out);

    return written_bytes;
}

uint64_t PreciseHashAnnotator::serialize(const std::string &filename) const {
    std::ofstream fout(filename);
    return serialize(fout);
}

void PreciseHashAnnotator::load(std::istream &in) {
    size_t num_prefix_cols = serialization::loadNumber(in);

    prefix_indices_.clear();
    while (num_prefix_cols--) {
        prefix_indices_.insert(serialization::loadNumber(in));
    }

    annotation_exact.load(in);
}

void PreciseHashAnnotator::load(const std::string &filename) {
    std::ifstream fin(filename);
    load(fin);
    fin.close();
}

uint64_t PreciseHashAnnotator::export_rows(std::ostream &out, bool permute) const {
    uint64_t written_bytes = serialization::serializeNumber(out, annotation_exact.kmer_map_.size());
    std::unordered_map<size_t, size_t> index_map;
    if (permute && prefix_indices_.size()) {
        index_map = compute_permutation_map();
    }
    /*
    for (auto &kmer : annotation_exact.kmer_map_) {
        auto annot = annotation_from_kmer(kmer.first);
        written_bytes += serialization::serializeNumberVector(out,
                index_map.size() ? permute_indices(annot, index_map) : annot);
    }
    */
    for (size_t i = 0; i < size(); ++i) {
        auto annot = annotate_edge(i, permute);
        written_bytes += serialization::serializeNumberVector(out, annot);
    }
    return written_bytes;
}

uint64_t PreciseHashAnnotator::export_rows(const std::string &filename, bool permute) const {
    std::ofstream fout(filename);
    return export_rows(fout, permute);
}

void PreciseHashAnnotator::add_sequence(const std::string &sequence,
                                        pos_t column,
                                        bool rooted) {
    std::string preprocessed_seq = graph_.transform_sequence(sequence, rooted);

    // Don't annotate short sequences
    if (preprocessed_seq.size() < graph_.get_k() + 1)
        return;

    if (column < static_cast<pos_t>(-1) && column >= annotation_exact.size())
        annotation_exact.resize(column + 1);

    for (size_t i = 0; i + graph_.get_k() < preprocessed_seq.size(); ++i) {
        annotation_exact.insert(
            &preprocessed_seq[i],
            &preprocessed_seq[i] + graph_.get_k() + 1, column
        );
    }
}

void PreciseHashAnnotator::add_column(const std::string &sequence, bool rooted) {
    add_sequence(sequence, annotation_exact.size(), rooted);
}

std::vector<uint64_t>
PreciseHashAnnotator::annotation_from_kmer(const std::string &kmer, bool permute) const {
    assert(kmer.length() == graph_.get_k() + 1);
    auto annot = annotation_exact.find(kmer.data(), kmer.data() + kmer.size());
    if (!permute || prefix_indices_.empty())
        return annot;

    return permute_indices(annot, compute_permutation_map());
}

size_t compute_optimal_num_hashes(double bloom_fpp, double bloom_size_factor = -1) {
    return bloom_size_factor > -0.5
           ? static_cast<size_t>(ceil(bloom_size_factor * std::log(2)))
           : (bloom_fpp > -0.5
                 ? static_cast<size_t>(ceil(-std::log2(bloom_fpp)))
                 : 0);
}

double compute_optimal_bloom_size_factor(double bloom_fpp) {
    return -std::log2(bloom_fpp) / std::log(2);
}

// Computes optimal `bloom_size_factor` and `num_hash_functions` automatically
BloomAnnotator::BloomAnnotator(const DeBruijnGraphWrapper &graph,
                               double bloom_fpp,
                               bool verbose)
      : graph_(graph),
        bloom_size_factor_(compute_optimal_bloom_size_factor(bloom_fpp)),
        bloom_fpp_(bloom_fpp),
        annotation(compute_optimal_num_hashes(bloom_fpp_, bloom_size_factor_)),
        total_traversed_(0),
        verbose_(verbose) {
    if (!annotation.num_hash_functions()) {
        std::cerr << "ERROR: invalid Bloom filter parameters" << std::endl;
        exit(1);
    }
}

// If not provided, computes optimal `num_hash_functions` automatically
BloomAnnotator::BloomAnnotator(const DeBruijnGraphWrapper &graph,
                               double bloom_size_factor,
                               size_t num_hash_functions,
                               bool verbose)
      : graph_(graph),
        bloom_size_factor_(bloom_size_factor),
        bloom_fpp_(exp(-bloom_size_factor_ * std::log(2) * std::log(2))),
        annotation(num_hash_functions
                      ? num_hash_functions
                      : compute_optimal_num_hashes(bloom_fpp_, bloom_size_factor_)),
        total_traversed_(0),
        verbose_(verbose) {
    if (!annotation.num_hash_functions()) {
        std::cerr << "ERROR: invalid Bloom filter parameters" << std::endl;
        exit(1);
    }
}

size_t BloomAnnotator::num_hash_functions() const {
    return annotation.num_hash_functions();
}

double BloomAnnotator::size_factor() const {
    return bloom_size_factor_;
}

double BloomAnnotator::approx_false_positive_rate() const {
    return bloom_fpp_;
}

size_t BloomAnnotator::get_size(size_t i) const {
    return annotation[i].size();
}

void BloomAnnotator::add_sequence(const std::string &sequence, size_t column, size_t num_elements) {
    std::string preprocessed_seq = graph_.encode_sequence(sequence);

    // Don't annotate short sequences
    if (preprocessed_seq.size() < graph_.get_k() + 1)
        return;

    if (column >= annotation.size())
        annotation.resize(column + 1);

    if (annotation[column].size() == 0) {
        annotation[column].resize(static_cast<size_t>(
            bloom_size_factor_
            * static_cast<double>(num_elements ? num_elements
                                               : preprocessed_seq.size() - graph_.get_k())
            + 1
        ));
        if (annotation[column].size() == 0) {
            std::cerr << "ERROR: resize failed" << std::endl;
            exit(1);
        }
    }

    for (auto hash_it = CyclicHashIterator(preprocessed_seq,
                                           graph_.get_k() + 1,
                                           annotation.num_hash_functions());
                !hash_it.is_end(); ++hash_it) {
        annotation.insert(*hash_it, column);
    }
}

void BloomAnnotator::add_column(const std::string &sequence, size_t num_elements) {
    add_sequence(sequence, annotation.size(), num_elements);
}

std::vector<uint64_t>
BloomAnnotator::annotation_from_kmer(const std::string &kmer) const {
    return annotation.find(annotation.compute_hash(kmer));
}

std::vector<uint64_t>
BloomAnnotator::get_annotation(DeBruijnGraphWrapper::edge_index i) const {
    return annotation_from_kmer(kmer_from_index(i));
}

std::vector<uint64_t>
BloomAnnotator::get_annotation_corrected(DeBruijnGraphWrapper::edge_index i,
                                         bool check_both_directions,
                                         size_t path_cutoff) const {
    //initial raw annotation
    std::string orig_kmer = kmer_from_index(i);
    assert(orig_kmer.length() == graph_.get_k() + 1);
    auto hasher = CyclicMultiHash(orig_kmer, annotation.num_hash_functions());

    //auto curannot = annotation_from_kmer(orig_kmer);
    auto curannot = annotation.find(hasher.get_hash());

    // Dummy edges are not supposed to be annotated
    if (graph_.is_dummy_edge(orig_kmer)) {
        curannot.assign(curannot.size(), 0);
        return curannot;
    }

    size_t pcount_old = hash_annotate::popcount(curannot);

    if (!pcount_old)
        return curannot;

    char cur_edge = orig_kmer.back();
    auto j = i;
    size_t path = 0;
    while (path++ < path_cutoff) {
        const_cast<size_t&>(total_traversed_)++;

        //traverse forward
        j = graph_.next_edge(j, cur_edge);

        cur_edge = graph_.get_edge_label(j);

        //check outdegree
        if (graph_.is_dummy_label(cur_edge)
                || !graph_.has_the_only_outgoing_edge(j)
                || (check_both_directions && !graph_.has_the_only_incoming_edge(j)))
            break;

        hasher.update(cur_edge);

        //bitwise AND annotations
        auto nextannot = hash_annotate::merge_and(
            curannot,
            annotation.find(hasher.get_hash())
        );

        //check popcounts
        size_t pcount_new = hash_annotate::popcount(nextannot);

        assert(pcount_new <= pcount_old);

        //if zero, then start of new sequence
        if (pcount_new == 0)
            break;

        //path length stopping conditions
        if (pcount_new < pcount_old) {
            curannot = nextannot;
            path = 0;
            pcount_old = pcount_new;
        }
    }

    //backward correction
    std::vector<DeBruijnGraphWrapper::edge_index> indices(graph_.get_k() + 1);
    assert(orig_kmer.length() == indices.size());
    indices[0] = i;

    auto back_hasher = CyclicMultiHash(orig_kmer, annotation.num_hash_functions());
    j = i;
    for (size_t m = 0; m < graph_.get_k(); ++m) {
        j = graph_.prev_edge(j);
        indices[m + 1] = j;
    }
    size_t back = graph_.get_k(); //index of the back
    assert(orig_kmer.front() == graph_.get_edge_label(indices[back]));
    assert(orig_kmer.back() == graph_.get_edge_label(indices[(back + 1) % indices.size()]));
    path = 0;
    while (graph_.has_the_only_incoming_edge(indices[(back + 1) % indices.size()])
            && (!check_both_directions
                || graph_.has_the_only_outgoing_edge(indices[(back + 1) % indices.size()]))
            && path++ < path_cutoff) {
        const_cast<size_t&>(total_traversed_)++;

        indices[(back + 1) % indices.size()] = graph_.prev_edge(indices[back]);
        back = (back + 1) % indices.size();

        char cur_first = graph_.get_edge_label(indices[back]);

        if (graph_.is_dummy_label(cur_first))
            break;

        back_hasher.reverse_update(cur_first);

        auto nextannot = hash_annotate::merge_and(
            curannot,
            annotation.find(back_hasher.get_hash())
        );

        auto pcount_new = hash_annotate::popcount(nextannot);

        assert(pcount_new <= pcount_old);

        if (pcount_new == 0)
            break;

        if (pcount_new < pcount_old) {
            curannot = nextannot;
            path = 0;
            pcount_old = pcount_new;
        }
    }

	return curannot;
}

void BloomAnnotator::test_fp_all(const PreciseAnnotator &annotation_exact,
                                 size_t num,
                                 bool check_both_directions) const {
    double fp_per_bit = 0;
    double fp_pre_per_bit = 0;
    double fn_per_bit = 0;
    uint64_t fp = 0;
    uint64_t fp_pre = 0;
    uint64_t fn = 0;
    uint64_t total = 0;
    assert(num);
    size_t step = std::max(
        static_cast<size_t>(1),
        static_cast<size_t>((graph_.last_edge() - graph_.first_edge() + 1) / num)
    );
    for (auto i = graph_.first_edge(); i <= graph_.last_edge(); i += step) {
        if (graph_.is_dummy_edge(kmer_from_index(i)))
            continue;
        total++;
        auto stats = test_fp(i, annotation_exact, check_both_directions);
        fp_pre_per_bit += (double)stats[0];
        fp_per_bit     += (double)stats[1];
        fn_per_bit     += (double)stats[2];
        fp_pre         += stats[0] > 0;
        fp             += stats[1] > 0;
        fn             += stats[2] > 0;
        //if (stats[2] > 0)
        //    step -= step - 1;
    }
    std::cout << "\n";
    std::cout << "Total:\t" << total << "\n";
    std::cout << "Post:\t"
              << "FP(edges):\t" << fp << "\t"
              << "FP(per edge):\t" << (double)fp / (double)total << "\t"
              << "FN(edges):\t" << fn
              << "\n";
    std::cout << "Pre:\t"
              << "FP(edges):\t" << fp_pre << "\t"
              << "FP(per edge):\t" << (double)fp_pre / (double)total << "\t"
              << "\n";
    std::cout << "Per bit" << "\n";
    std::cout << "Post:\t"
              << "FP(bits/edge):\t" << (double)fp_per_bit / (double)total << "\t"
              << "Avg. FPP:\t" << (double)fp_per_bit / (double)total / (double)annotation.size() << "\t"
              << "FN(bits):\t" << (double)fn_per_bit
              << "\n";
    std::cout << "Pre:\t"
              << "FP(bits/edge):\t" << (double)fp_pre_per_bit / (double)total << "\t"
              << "Avg. FPP:\t" << (double)fp_pre_per_bit / (double)total / (double)annotation.size() << "\t"
              << "\n";
    std::cout << "Total traversed: " << total_traversed_ << "\n";
}

uint64_t BloomAnnotator::serialize(std::ostream &out) const {
    return annotation.serialize(out);
}

uint64_t BloomAnnotator::serialize(const std::string &filename) const {
    std::ofstream out(filename);
    return serialize(out);
}

bool BloomAnnotator::load(std::istream &in) {
    try {
        annotation.load(in);
    } catch (...) {
        return false;
    }
    return true;
}

bool BloomAnnotator::load(const std::string &filename) {
    std::ifstream in(filename);
    if (!in.good() || !load(in))
        return false;
    in.close();
    return true;
}

std::vector<size_t>
BloomAnnotator::unpack(const std::vector<uint64_t> &packed) {
    std::vector<size_t> labels;
    for (size_t i = 0; i < packed.size() * 64; ++i) {
        if (packed[i / 64] & (1llu << (i % 64)))
            labels.push_back(i);
    }
    return labels;
}

std::string
BloomAnnotator::kmer_from_index(DeBruijnGraphWrapper::edge_index index) const {
    assert(index <= graph_.last_edge());
    assert(index >= graph_.first_edge());
    return graph_.get_node_kmer(index) + graph_.get_edge_label(index);
}

std::string
PreciseHashAnnotator::get_kmer(DeBruijnGraphWrapper::edge_index i) const {
    return graph_.get_node_kmer(i) + graph_.get_edge_label(i);
}

std::vector<uint64_t>
PreciseHashAnnotator::annotate_edge(DeBruijnGraphWrapper::edge_index i, bool permute) const {
    return annotation_from_kmer(get_kmer(i), permute);
}

std::set<pos_t>
PreciseHashAnnotator::annotate_edge_indices(DeBruijnGraphWrapper::edge_index i, bool permute) const {
    auto find = annotation_exact.kmer_map_.find(get_kmer(i));
    if (find == annotation_exact.kmer_map_.end())
        return {};
    if (!permute)
        return find->second;
    auto index_map = compute_permutation_map();
    if (index_map.empty())
        return find->second;
    std::vector<pos_t> find_mapped_v;
    find_mapped_v.reserve(find->second.size());
    std::transform(find->second.begin(),
                   find->second.end(),
                   std::back_inserter(find_mapped_v),
                   [&](pos_t i) {
                       return index_map[i];
                   });
    return std::set<pos_t>(find_mapped_v.begin(), find_mapped_v.end());
    /*
    std::set<pos_t> find_mapped;
    std::transform(find->second.begin(),
                   find->second.end(),
                   std::inserter(find_mapped, find_mapped.begin()),
                   [&](pos_t i) {
                       return index_map[i];
                   });
    return find_mapped;
    */
}

std::vector<uint64_t>
BloomAnnotator::test_fp(DeBruijnGraphWrapper::edge_index i,
                        const PreciseAnnotator &annotation_exact,
                        bool check_both_directions) const {

    auto int_kmer = kmer_from_index(i);

    auto test = annotation_from_kmer(int_kmer);

    auto test_exact = annotation_exact.annotate_edge(i);

    auto curannot = get_annotation_corrected(i, check_both_directions);

    auto jt = test.begin();
    auto kt = test_exact.begin();
    auto lt = curannot.begin();

    std::vector<uint64_t> stats(3, 0);

    for (; jt != test.end(); ++jt, ++kt, ++lt) {
        //check for false negatives
        if ((*jt | *kt) != *jt) {
            std::cerr << "Encoding " << i << " failed" << std::endl;
            std::cout << "FN: " << int_kmer << std::endl;
            auto unpacked_labels = unpack(test_exact);
            std::cout << "True annotation:\t";
            for (size_t value : unpacked_labels) {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            unpacked_labels = unpack(test);
            std::cout << "Uncorrected annotation:\t";
            for (size_t value : unpacked_labels) {
                std::cout << value << " ";
            }
            std::cout << std::endl;

            //exit(1);
        }
        //correction introduced extra bits
        if ((*lt | *jt) != *jt) {
            std::cerr << "False bits added" << std::endl;
            //exit(1);
        }
        //false positives before correction
        if ((*jt | *kt) != *kt) {
            //stats[0] += __builtin_popcountll(*jt) - __builtin_popcountll(*kt);
            stats[0] += static_cast<uint64_t>(__builtin_popcountll(*jt & (~(*kt))));
        }
        //false positives after correction
        if ((*lt | *kt) != *kt) {
            //stats[1] += __builtin_popcountll(*lt) - __builtin_popcountll(*kt);
            stats[1] += static_cast<uint64_t>(__builtin_popcountll(*lt & (~(*kt))));
            if (verbose_)
                std::cout << "FP: " << int_kmer << std::endl;
        }
        //false negatives after correction
        if ((*lt | *kt) != *lt) {
            //stats[2] += __builtin_popcountll(*lt | *kt) - __builtin_popcountll(*lt);
            stats[2] += static_cast<uint64_t>(__builtin_popcountll(*kt & (~(*lt))));
            if (verbose_) {
                std::cout << "FN: " << int_kmer << std::endl;
                auto unpacked_labels = unpack(test_exact);
                std::cout << "True annotation:\t";
                for (size_t value : unpacked_labels) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
                unpacked_labels = unpack(curannot);
                std::cout << "Corrected annotation:\t";
                for (size_t value : unpacked_labels) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
            }
        }
        //if (stats[0] && stats[1] && stats[2]) {
        //    break;
        //}
    }
    return stats;
}

} // namespace hash_annotate
