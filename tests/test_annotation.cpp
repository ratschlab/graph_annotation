#include <random>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include "gtest/gtest.h"
#include "dbg_bloom_annotator.hpp"
#include "wavelet_trie_annotator.hpp"
#include "serialization.hpp"
#include "hashers.hpp"
#include "dbg_hash.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const size_t num_random_kmers = 1000;

struct uint256_t {
    uint64_t m_lo;
    uint64_t m_mid;
    __uint128_t m_high;
};

const std::vector<size_t> num_threads = {1, 4};

std::vector<std::string> generate_kmers(size_t num, size_t k = 82) {
    std::vector<std::string> kmers(num);
    for (auto &kmer : kmers) {
        for (size_t j = 0; j < k; ++j) {
            kmer.push_back((rand() % 26) + 65);
        }
    }
    return kmers;
}

TEST(Annotate, RandomTestNoFalseNegative) {
    //create annotation
    hash_annotate::BloomHashAnnotation bloom_anno(7);
    hash_annotate::BloomFilter bloom(num_random_kmers);
    hash_annotate::ExactHashAnnotation exact_anno;
    hash_annotate::ExactFilter exact;
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    size_t total = 0, fp = 0;
    for (size_t i = 0; i < kmers.size(); ++i) {
        if (i < kmers.size() / 2) {
            bloom.insert(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()));
            exact.insert(exact_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()));
            ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
            ASSERT_TRUE(exact.find(exact_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
        } else {
            if (exact.find(exact_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()))) {
                ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
            }
            if (!exact.find(exact_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()))) {
                total++;
                if (bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()))) {
                    fp++;
                }
            }
        }
    }
    EXPECT_TRUE(total >= fp) << "Total: " << total << " FP: " << fp << std::endl;
}

TEST(Annotate, BloomFilterEmpty) {
    hash_annotate::BloomHashAnnotation bloom_anno(7);
    hash_annotate::BloomFilter empty(num_random_kmers);

    EXPECT_EQ(0, empty.occupancy());
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        ASSERT_FALSE(empty.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
    }
}

TEST(Annotate, BloomFilterEmptyHash) {
    hash_annotate::BloomHashAnnotation bloom_anno;
    hash_annotate::BloomFilter bloom(num_random_kmers);
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom.insert(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()));
        ASSERT_FALSE(bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
    }
}

TEST(Annotate, BloomEquality) {
    hash_annotate::BloomHashAnnotation bloom_anno(7);
    hash_annotate::BloomFilter bloom(num_random_kmers);
    hash_annotate::BloomFilter empty1(1), empty2(num_random_kmers);
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom.insert(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()));
        ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
        ASSERT_FALSE(empty1.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
        ASSERT_FALSE(empty2.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
    }
    EXPECT_EQ(bloom, bloom);
    EXPECT_NE(empty1, bloom);
    EXPECT_NE(empty2, bloom);
}

TEST(Annotate, Serialize) {
    hash_annotate::BloomHashAnnotation bloom_anno(7);
    hash_annotate::BloomFilter bloom(num_random_kmers);
    hash_annotate::ExactHashAnnotation exact_anno;
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom.insert(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length()));
        ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(kmers[i].data(), kmers[i].data() + kmers[i].length())));
        exact_anno.insert(kmers[i].data(), kmers[i].data() + kmers[i].length(), 0);
    }

    std::ofstream outstream(test_dump_basename + "_bloomser");
    bloom.serialize(outstream);
    outstream.close();

    outstream.open(test_dump_basename + "_exactser");
    exact_anno.serialize(outstream);
    outstream.close();

    std::ifstream instream(test_dump_basename + "_bloomser");
    hash_annotate::BloomFilter bloom_alt(1);
    bloom_alt.load(instream);
    EXPECT_EQ(bloom, bloom_alt);
    instream.close();

    instream.open(test_dump_basename + "_exactser");
    hash_annotate::ExactHashAnnotation exact_alt;
    exact_alt.load(instream);
    EXPECT_EQ(exact_anno, exact_alt);
}

TEST(Annotate, BloomHashSerialize) {
    hash_annotate::BloomHashAnnotation bloom_anno(7);
    bloom_anno.append_bit(num_random_kmers);
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom_anno.insert(kmers[i].data(), kmers[i].data() + kmers[i].length(), 0);
        ASSERT_EQ(1u, bloom_anno.find(kmers[i].data(), kmers[i].data() + kmers[i].length(), 0)[0]);
    }
    std::ofstream outstream(test_dump_basename + "_bloomser");
    bloom_anno.serialize(outstream);
    outstream.close();

    std::ifstream instream(test_dump_basename + "_bloomser");
    hash_annotate::BloomHashAnnotation bloom_alt;
    bloom_alt.load(instream);
    EXPECT_EQ(bloom_anno, bloom_alt);
}


TEST(Annotate, RandomHashAnnotator) {
    hash_annotate::BloomHashAnnotation bloomhash(7);
    hash_annotate::ExactHashAnnotation exacthash;
    size_t num_bits = 5;
    std::vector<size_t> bounds(num_bits);
    std::iota(bounds.begin(), bounds.end(), 0);
    for (size_t i = 0; i < num_bits; ++i) {
        bloomhash.append_bit(1000);
        exacthash.append_bit();
    }
    ASSERT_EQ(bloomhash.size(), num_bits);
    ASSERT_EQ(exacthash.size(), num_bits);
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        size_t pick_bits = 0;
        if (i < kmers.size()) {
            //insert into random bit positions
            pick_bits = rand() % (1u << num_bits);
            ASSERT_TRUE(pick_bits < (1u << num_bits));
            for (size_t j = 0; j < num_bits; ++j) {
                if (((1lu << j) | pick_bits) != pick_bits)
                    continue;

                //insert
                auto testbloom = bloomhash.insert(kmers[i].data(), kmers[i].data() + kmers[i].length(), j);
                auto testexact = exacthash.insert(kmers[i].data(), kmers[i].data() + kmers[i].length(), j);
                //check if it's there
                testbloom = bloomhash.find(kmers[i].data(), kmers[i].data() + kmers[i].length(), j);
                testexact = exacthash.find(kmers[i].data(), kmers[i].data() + kmers[i].length(), j);

                //test OR
                auto testbloom_merged = testbloom;
                hash_annotate::merge_or(testbloom_merged, testexact);
                ASSERT_TRUE(hash_annotate::equal(testbloom, testbloom_merged));

                //test bit
                ASSERT_TRUE(hash_annotate::test_bit(testbloom, j));
                ASSERT_TRUE(hash_annotate::test_bit(testexact, j));

                //test AND
                auto testbloom_and = testbloom;
                hash_annotate::merge_and(testbloom_merged, testexact);
                ASSERT_TRUE(hash_annotate::equal(testexact, testbloom_and));
            }
        }
        auto testbloom = bloomhash.find(kmers[i].data(), kmers[i].data() + kmers[i].length());
        auto testexact = exacthash.find(kmers[i].data(), kmers[i].data() + kmers[i].length());
        ASSERT_TRUE(hash_annotate::equal(testbloom, hash_annotate::merge_or(testbloom, testexact)));
    }
}

TEST(Annotate, HashIterator) {
    std::string test_string;
    for (size_t i = 0; i < 8; ++i) {
        test_string += std::string("$NATGC");
    }
    ASSERT_EQ(48llu, test_string.length());

    size_t num_hash_functions = 5;
    size_t kmer_size = 20;

    hash_annotate::CyclicHashIterator hash_it(test_string, kmer_size, num_hash_functions);
    hash_annotate::CyclicMultiHash hash_it_up(test_string.substr(0, kmer_size),
                                              num_hash_functions);
    ASSERT_EQ(num_hash_functions, hash_it->size());

    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        ASSERT_FALSE(hash_it.is_end());
        ASSERT_NE('\0', test_string[i + kmer_size - 1]);
        for (uint32_t j = 0; j < num_hash_functions; ++j) {
            ASSERT_EQ(hash_it->at(j), hash_it_up.get_hash()[j]);
        }
        ++hash_it;
        hash_it_up.update(test_string[i + kmer_size]);
    }
    ASSERT_TRUE(hash_it.is_end());
}

TEST(Annotate, HashIteratorEmpty) {
    size_t num_hash_functions = 1;
    hash_annotate::CyclicHashIterator hash_it("", 1, num_hash_functions);
    ASSERT_TRUE(hash_it.is_end());

    hash_annotate::CyclicHashIterator hash_it2("1", 5, num_hash_functions);
    ASSERT_TRUE(hash_it.is_end());
}

TEST(Annotate, BloomFilter) {
    for (size_t k = 10; k <= 90; k += 10) {
        auto kmers = generate_kmers(num_random_kmers, k + 1);
        size_t num_seqs = 10;
        size_t size_chunk = kmers.size() / num_seqs;

        std::vector<std::string> sequences(num_seqs);
        for (size_t i = 0; i < num_seqs; ++i) {
            sequences[i] = std::accumulate(
                    kmers.begin() + i * size_chunk,
                    kmers.begin() + (i + 1) * size_chunk,
                    std::string(""));
        }

        DBGHash graph(k);
        hash_annotate::PreciseHashAnnotator precise(graph);

        for (size_t i = 0; i < num_seqs; ++i) {
            graph.add_sequence(sequences[i]);
            precise.add_sequence(sequences[i], i);
        }

        ASSERT_EQ(graph.get_num_edges(), precise.size());

        for (double bloom_fpp = 0.05; bloom_fpp < 1.0; bloom_fpp *= 3) {
            hash_annotate::BloomAnnotator bloom(graph, bloom_fpp);
            for (size_t i = 0; i < num_seqs; ++i) {
                bloom.add_sequence(sequences[i], i);
            }
            for (auto &kmer : kmers) {
                auto bloom_annot = bloom.annotation_from_kmer(kmer);
                auto precise_annot = precise.annotation_from_kmer(kmer);
                EXPECT_TRUE(hash_annotate::equal(
                            hash_annotate::merge_or(bloom_annot, precise_annot),
                            bloom_annot));
            }
        }
    }
}

TEST(Annotate, WaveletTrie) {
    for (size_t k = 10; k <= 90; k += 10) {
        auto kmers = generate_kmers(num_random_kmers, k + 1);
        size_t num_seqs = 10;
        size_t size_chunk = kmers.size() / num_seqs;

        std::vector<std::string> sequences(num_seqs);
        for (size_t i = 0; i < num_seqs; ++i) {
            sequences[i] = std::accumulate(
                    kmers.begin() + i * size_chunk,
                    kmers.begin() + (i + 1) * size_chunk,
                    std::string(""));
        }

        DBGHash graph(k);
        hash_annotate::PreciseHashAnnotator precise(graph);
        annotate::WaveletTrieAnnotator wts_ext(graph, 2);

        for (size_t i = 0; i < num_seqs; ++i) {
            graph.add_sequence(sequences[i]);
            precise.add_sequence(sequences[i], i);
            wts_ext.add_sequence(sequences[i], i);
            annotate::WaveletTrieAnnotator wts_pre(precise, graph);
            if (wts_pre != wts_ext) {
                wts_pre.print();
                std::cout << "\n";
                wts_ext.print();
                std::cout << "\n";
            }
            ASSERT_EQ(wts_pre, wts_ext) << i;
            /*
            for (size_t j = 0; j < wts_ext.size(); ++j) {
                if (!hash_annotate::equal(precise.annotate_edge(j, true),wts_ext.annotate_edge(j, true))) {
                    wts_ext.print();
                    std::cout << "\n";
                    annotate::WaveletTrieAnnotator wts_pr(precise, graph);
                    wts_pr.print();
                    std::cout << "\n";
                }
                ASSERT_TRUE(hash_annotate::equal(
                            precise.annotate_edge(j, true),
                            wts_ext.annotate_edge(j, true))) << i << "\t" << j;
            }
            */
        }

        ASSERT_EQ(graph.get_num_edges(), precise.size());
        ASSERT_EQ(graph.get_num_edges(), wts_ext.size());

        std::vector<annotate::WaveletTrieAnnotator> wtrs{std::move(wts_ext)};

        for (auto p : num_threads) {
            //std::cout << "Checking normal" << std::endl;
            annotate::WaveletTrieAnnotator wtr(precise, graph, p);
            EXPECT_EQ(graph.get_num_edges(), wtr.size());
            for (size_t i = 0; i < wtr.size(); ++i) {
                ASSERT_TRUE(hash_annotate::equal(
                            precise.annotate_edge(i, true),
                            wtr.annotate_edge(i))) << p << "\t" << i;
                auto kmer = graph.get_node_kmer(i) + graph.get_edge_label(i);
                ASSERT_TRUE(hash_annotate::equal(
                            precise.annotation_from_kmer(kmer, true),
                            wtr.annotation_from_kmer(kmer))) << p << "\t" << i;
            }
            wtrs.push_back(wtr);
            //std::cout << "Serializing" << std::endl;
            std::ofstream out(test_dump_basename + "_wtr");
            wtr.serialize(out);
            out.close();

            //std::cout << "Loading" << std::endl;
            std::ifstream in(test_dump_basename + "_wtr");
            annotate::WaveletTrieAnnotator wtr_out(graph, p);
            wtr_out.load(in);
            wtrs.push_back(wtr_out);
            in.close();


            //std::cout << "Dumping precise" << std::endl;
            out.open(test_dump_basename + "_precise");
            precise.serialize(out);
            out.close();

            in.open(test_dump_basename + "_precise");
            hash_annotate::PreciseHashAnnotator precise_file(graph);

            //std::cout << "Loading precise" << std::endl;
            precise_file.load(in);

            //std::cout << "Checking precise" << std::endl;
            ASSERT_EQ(precise, precise_file);
            in.close();

            //std::cout << "Construct from object" << std::endl;
            annotate::WaveletTrieAnnotator wtr_file(precise_file, graph, p);
            EXPECT_EQ(graph.get_num_edges(), wtr_file.size());
            for (size_t i = 0; i < wtr_file.size(); ++i) {
                ASSERT_TRUE(hash_annotate::equal(
                            precise_file.annotate_edge(i),
                            wtr_file.annotate_edge(i)));
            }
            wtrs.push_back(wtr_file);

            //std::cout << "Construct from file" << std::endl;
            in.open(test_dump_basename + "_precise");
            annotate::WaveletTrieAnnotator wtr_file2(graph, p);
            wtr_file2.load_from_precise_file(in, p);
            EXPECT_EQ(graph.get_num_edges(), wtr_file2.size());
            for (size_t i = 0; i < wtr_file2.size(); ++i) {
                ASSERT_TRUE(hash_annotate::equal(
                            precise_file.annotate_edge(i),
                            wtr_file2.annotate_edge(i)))
                    << "FAILED ON " << p << " threads\n"
                    << "INDEX " << i << "\n"
                    << precise_file.annotate_edge(i)[0] << "\t"
                    << wtr_file2.annotate_edge(i)[0] << "\t"
                    << "NumPrefix " << precise_file.num_prefix_columns()
                    << std::endl;
            }
            wtrs.push_back(wtr_file2);
            in.close();

            //std::cout << "Serialization" << std::endl;
            out.open(test_dump_basename + "_wtr");
            wtr_file.serialize(out);
            wtr_file2.serialize(out);
            out.close();
            in.open(test_dump_basename + "_wtr");
            annotate::WaveletTrieAnnotator wtr_file3(graph, p);
            annotate::WaveletTrieAnnotator wtr_file4(graph, p);
            wtr_file3.load(in);
            wtr_file4.load(in);
            in.close();
            wtrs.push_back(wtr_file3);
            wtrs.push_back(wtr_file4);
        }

        //std::cout << "Checking equality" << std::endl;
        for (size_t i = 1; i < wtrs.size(); ++i) {
            if (wtrs[i - 1] != wtrs[i]) {
                wtrs[i - 1].print();
                std::cout << "\n";
                wtrs[i].print();
                std::cout << "\n";
            }
            ASSERT_EQ(wtrs[i - 1], wtrs[i]) << i;
        }
    }
}

TEST(Annotate, ExportColsWithWithoutRearrange) {
    for (size_t k = 10; k < 90; k += 10) {
        auto kmers = generate_kmers(num_random_kmers, k + 1);
        size_t num_seqs = 10;
        size_t size_chunk = kmers.size() / num_seqs;
        std::vector<std::string> sequences(num_seqs);
        for (size_t i = 0; i < num_seqs; ++i) {
            sequences[i] = std::accumulate(
                    kmers.begin() + i * size_chunk,
                    kmers.begin() + (i + 1) * size_chunk,
                    std::string(""));
            EXPECT_LT(0u, sequences[i].length());
        }

        DBGHash graph(k);
        hash_annotate::PreciseHashAnnotator precise(graph);
        for (size_t i = 0; i < num_seqs; ++i) {
            graph.add_sequence(sequences[i]);
            precise.add_sequence(sequences[i], i);
        }
        EXPECT_EQ(graph.get_num_edges(), precise.size());
        annotate::WaveletTrieAnnotator wtr(precise, graph);
        EXPECT_EQ(precise.size(), wtr.size());

        // compute column permutation map
        std::vector<size_t> prefix_cols = {2, 5};
        precise.make_columns_prefix(prefix_cols.begin(), prefix_cols.end());
        auto index_map = precise.compute_permutation_map();
        ASSERT_GT(precise.num_columns(), *prefix_cols.rbegin());
        ASSERT_EQ(prefix_cols.size(), precise.num_prefix_columns());
        annotate::WaveletTrieAnnotator wtr_rearr(precise, graph);
        ASSERT_EQ(wtr.size(), wtr_rearr.size());

        // dump data
        std::ofstream os(test_dump_basename + "_precise");
        std::ofstream os_rearrange(test_dump_basename + "_rearrange");
        std::ofstream os_precise(test_dump_basename + "_precise2");
        std::ofstream os_wtr(test_dump_basename + "_wtr");
        std::ofstream os_wtrr(test_dump_basename + "_wtrr");
        precise.export_rows(os, false);
        precise.export_rows(os_rearrange);
        precise.serialize(os_precise);
        wtr.serialize(os_wtr);
        wtr_rearr.serialize(os_wtrr);
        os.close();
        os_rearrange.close();
        os_precise.close();
        os_wtr.close();
        os_wtrr.close();
        
        // reload data
        std::ifstream is(test_dump_basename + "_precise");
        std::ifstream is_rearrange(test_dump_basename + "_rearrange");
        std::ifstream is_precise(test_dump_basename + "_precise2");
        std::ifstream is_wtr(test_dump_basename + "_wtr");
        std::ifstream is_wtrr(test_dump_basename + "_wtrr");
        hash_annotate::PreciseHashAnnotator precise_file(graph);
        precise_file.load(is_precise);
        ASSERT_EQ(precise, precise_file);
        annotate::WaveletTrieAnnotator wtr_file(graph);
        wtr_file.load(is_wtr);
        annotate::WaveletTrieAnnotator wtr_filer(graph);
        wtr_filer.load(is_wtrr);
        is_precise.close();
        ASSERT_EQ(precise.num_prefix_columns(), precise_file.num_prefix_columns());
        ASSERT_EQ(precise.num_columns(), precise_file.num_columns());
        ASSERT_EQ(wtr.size(), wtr_file.size());
        ASSERT_EQ(wtr.size(), wtr_filer.size());

        // check data
        size_t num_kmers = serialization::loadNumber(is);
        size_t num_kmers_rearrange = serialization::loadNumber(is_rearrange);
        ASSERT_EQ(graph.get_num_edges(), num_kmers);
        ASSERT_EQ(graph.get_num_edges(), num_kmers_rearrange);
        size_t wtr_i = 0;
        ASSERT_EQ(precise.size(), num_kmers);
        for (size_t j = 0; j < precise.size(); ++j) {
            auto cur_annot = serialization::loadNumberVector(is);
            auto cur_annot_rearrange = serialization::loadNumberVector(is_rearrange);
            auto ref_annot = precise.annotate_edge(j, true);
            //auto ref_annot = precise.annotation_from_kmer(kmer.first, true);
            auto ref_annot_file = precise_file.annotate_edge(j, true);
            //auto ref_annot_file = precise_file.annotation_from_kmer(kmer.first, true);
            auto wtr_filer_annot = wtr_filer.annotate_edge(wtr_i);
            ASSERT_TRUE(hash_annotate::equal(ref_annot, cur_annot_rearrange));
            ASSERT_TRUE(hash_annotate::equal(ref_annot, ref_annot_file));
            ASSERT_TRUE(hash_annotate::equal(wtr.annotate_edge(wtr_i), wtr_file.annotate_edge(wtr_i)));
            ASSERT_TRUE(hash_annotate::equal(cur_annot, wtr_file.annotate_edge(wtr_i)));
            ASSERT_TRUE(hash_annotate::equal(ref_annot, wtr_filer.annotate_edge(wtr_i)));
            wtr_i++;
            for (auto &curi : index_map) {
                ASSERT_EQ(hash_annotate::test_bit(ref_annot, curi.second),
                          hash_annotate::test_bit(cur_annot, curi.first));
            }

            //check if each row is present in at least one category
            size_t i = 0;
            for (auto limb : ref_annot) {
                if (limb)
                    break;
                i++;
            }
            ASSERT_GT(ref_annot.size(), i);
        }
    }
}

//TODO: WRITE A UNIT TEST TO MAKE SURE BLOOM FILTERS ARE SUPERSET OF EXACT FILTER
/*
TEST(Annotate, HashIteratorInsert) {
    std::string test_string;
    for (size_t i = 0; i < 8; ++i) {
        test_string += std::string("$NATGC");
    }
    ASSERT_EQ(48llu, test_string.length());

    size_t num_hash_functions = 5;
    size_t kmer_size = 20;

    hash_annotate::MurmurHashIterator hash_it(test_string, num_hash_functions, kmer_size);
    ntHashIterator hash_nt_it(test_string, num_hash_functions, kmer_size);
    ASSERT_EQ(num_hash_functions, hash_it.size());

    hash_annotate::HashAnnotation<hash_annotate::BloomFilter> bloomhash(num_hash_functions);
    hash_annotate::HashAnnotation<hash_annotate::ExactFilter> exacthash;
    hash_annotate::HashAnnotation<hash_annotate::BloomFilter> bloomhash_it(num_hash_functions);
    hash_annotate::HashAnnotation<hash_annotate::ExactFilter> exacthash_it;
    hash_annotate::HashAnnotation<hash_annotate::ExactFilter> exacthash_it2;
    hash_annotate::HashAnnotation<hash_annotate::ExactFilter> exacthash_nt_it;

    bloomhash.append_bit(1000);
    exacthash.append_bit(1000);
    bloomhash_it.append_bit(1000);
    exacthash_it.append_bit(1000);
    exacthash_it2.append_bit(1000);
    exacthash_nt_it.append_bit(1000);

    while (hash_it != hash_it.end()) {
        exacthash_it2.insert(hash_annotate::MultiHash(*hash_it, num_hash_functions), 0);
        ++hash_it;
    }

    while (hash_nt_it != hash_nt_it.end()) {
        exacthash_nt_it.insert(hash_annotate::MultiHash(*hash_nt_it, num_hash_functions), 0);
        ++hash_nt_it;
    }

    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        bloomhash.insert(&test_string[i], &test_string[i] + kmer_size, 0);
        exacthash.insert(&test_string[i], &test_string[i] + kmer_size, 0);
        bloomhash_it.insert(hashes[i], 0);
        exacthash_it.insert(hashes[i], 0);
    }
    EXPECT_TRUE(bloomhash == bloomhash_it);
    EXPECT_TRUE(exacthash == exacthash_it);
    EXPECT_TRUE(exacthash == exacthash_it2);

    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        auto hash = hash_annotate::hash_murmur(std::string(&test_string[i], kmer_size), num_hash_functions, kmer_size)[0];
        ASSERT_TRUE(exacthash.find(hash)[0]);
        ASSERT_TRUE(exacthash_it.find(hash)[0]);
        ASSERT_TRUE(exacthash_it2.find(hash)[0]);
        ASSERT_TRUE(bloomhash.find(hash)[0]);
        ASSERT_TRUE(bloomhash_it.find(hash)[0]);
        //ASSERT_TRUE(bloomhash_it2.find(hash)[0]);
    }
}
*/
