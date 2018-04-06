#include <random>
#include <fstream>

#include "gtest/gtest.h"
#include "dbg_bloom_annotator.hpp"
#include "hashers.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const size_t num_random_kmers = 1000;

struct uint256_t {
    uint64_t m_lo;
    uint64_t m_mid;
    __uint128_t m_high;
};

std::vector<uint256_t> generate_kmers(size_t num) {
    std::vector<uint256_t> kmers(num);
    for (size_t i = 0; i < kmers.size(); ++i) {
        *(reinterpret_cast<uint64_t*>(&kmers[i]))     = rand();
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 1) = rand();
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 2) = rand();
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 3) = rand();
    }
    return kmers;
}

TEST(Annotate, RandomTestNoFalseNegative) {
    //create annotation
    hash_annotate::BloomHashAnnotation bloom_anno;
    hash_annotate::ExactHashAnnotation exact_anno;
    hash_annotate::BloomFilter bloom(num_random_kmers);
    hash_annotate::ExactFilter exact;
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

TEST(Annotate, BloomSerialize) {
    hash_annotate::BloomHashAnnotation bloom_anno;
    hash_annotate::BloomFilter bloom(num_random_kmers);
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom.insert(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1));
        ASSERT_TRUE(bloom.find(bloom_anno.compute_hash(&kmers[i], &kmers[i] + 1)));
    }
    std::ofstream outstream(test_dump_basename + "_bloomser");
    bloom.serialize(outstream);
    outstream.close();

    std::ifstream instream(test_dump_basename + "_bloomser");
    hash_annotate::BloomFilter bloom_alt(1);
    bloom_alt.load(instream);
    EXPECT_EQ(bloom, bloom_alt);
}

TEST(Annotate, BloomHashSerialize) {
    hash_annotate::BloomHashAnnotation bloom_anno;
    bloom_anno.append_bit(num_random_kmers);
    //generate a bunch of kmers
    auto kmers = generate_kmers(num_random_kmers);
    for (size_t i = 0; i < kmers.size(); ++i) {
        bloom_anno.insert(&kmers[i], &kmers[i] + 1, 0);
        ASSERT_EQ(1, bloom_anno.find(&kmers[i], &kmers[i] + 1, 0)[0]);
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
                auto testbloom = bloomhash.insert(&kmers[i], &kmers[i] + 1, j);
                auto testexact = exacthash.insert(&kmers[i], &kmers[i] + 1, j);
                //check if it's there
                testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1, j);
                testexact = exacthash.find(&kmers[i], &kmers[i] + 1, j);

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
        auto testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1);
        auto testexact = exacthash.find(&kmers[i], &kmers[i] + 1);
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
