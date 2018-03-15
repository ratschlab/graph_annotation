#include <random>
#include <sdsl/wavelet_trees.hpp>

#include "gtest/gtest.h"
#include "dbg_bloom_annotator.hpp"
#include "hashers.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const size_t num_random_kmers = 1000;


std::vector<sdsl::uint256_t> generate_kmers(size_t num) {
    std::vector<sdsl::uint256_t> kmers(num);
    int mod = pow(num, 0.25);
    for (size_t i = 0; i < kmers.size(); ++i) {
        *(reinterpret_cast<uint64_t*>(&kmers[i]))     = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 1) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 2) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 3) = rand() % mod;
    }
    return kmers;
}

TEST(Annotate, RandomTestNoFalseNegative) {
    //create annotation
    annotate::BloomHashAnnotation bloom_anno;
    annotate::ExactHashAnnotation exact_anno;
    annotate::BloomFilter bloom(num_random_kmers);
    annotate::ExactFilter exact;
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    size_t total = 0, fp = 0;
    for (size_t i = 0; i < kmers.size(); ++i) {
        if (i < kmers.size() / 2) {
            bloom.insert(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1));
            exact.insert(exact_anno.compute_hash(&kmers[i], &kmers[i] + 1));
            ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1)));
            ASSERT_TRUE(exact.find(exact_anno.compute_hash(&kmers[i], &kmers[i] + 1)));
        } else {
            if (exact.find(exact_anno.compute_hash(&kmers[i], &kmers[i] + 1))) {
                ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1)));
            }
            if (!exact.find(exact_anno.compute_hash(&kmers[i], &kmers[i] + 1))) {
                total++;
                if (bloom.find(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1))) {
                    fp++;
                }
            }
        }
    }
    EXPECT_TRUE(total >= fp) << "Total: " << total << " FP: " << fp << std::endl;
}

TEST(Annotate, RandomHashAnnotator) {
    annotate::BloomHashAnnotation bloomhash(7);
    annotate::ExactHashAnnotation exacthash;
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
                auto testbloom = bloomhash.insert(&kmers[i], &kmers[i] + 1, j);
                auto testexact = exacthash.insert(&kmers[i], &kmers[i] + 1, j);
                //check if it's there
                testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1, j);
                testexact = exacthash.find(&kmers[i], &kmers[i] + 1, j);

                //test OR
                auto testbloom_merged = testbloom;
                annotate::merge_or(testbloom_merged, testexact);
                ASSERT_TRUE(annotate::equal(testbloom, testbloom_merged));

                //test bit
                ASSERT_TRUE(annotate::test_bit(testbloom, j));
                ASSERT_TRUE(annotate::test_bit(testexact, j));

                //test AND
                auto testbloom_and = testbloom;
                annotate::merge_and(testbloom_merged, testexact);
                ASSERT_TRUE(annotate::equal(testexact, testbloom_and));
            }
        }
        auto testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1);
        auto testexact = exacthash.find(&kmers[i], &kmers[i] + 1);
        ASSERT_TRUE(annotate::equal(testbloom, annotate::merge_or(testbloom, testexact)));
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

    annotate::CyclicHashIterator hash_it(test_string, kmer_size, num_hash_functions);
    annotate::CyclicMultiHash hash_it_up(test_string.substr(0, kmer_size),
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
    annotate::CyclicHashIterator hash_it("", 1, num_hash_functions);
    ASSERT_TRUE(hash_it.is_end());

    annotate::CyclicHashIterator hash_it2("1", 5, num_hash_functions);
    ASSERT_TRUE(hash_it.is_end());
}

