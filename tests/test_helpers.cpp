#include <fstream>

#include "gtest/gtest.h"
#include "wavelet_trie.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

annotate::cpp_int generate_number() {
    annotate::cpp_int l_int = -1llu - 1;
    l_int <<= 64;
    l_int += -1llu - 1;
    return l_int;
}

TEST(CPPINT, Serialization) {
    auto l_int = generate_number();
    std::ofstream out(test_dump_basename + ".gmp");
    ASSERT_TRUE(out.good());
    annotate::serialize(out, l_int);
    out.close();

    std::ifstream in(test_dump_basename + ".gmp");
    ASSERT_TRUE(out.good());
    auto l_int_file = annotate::load(in);
    EXPECT_EQ(l_int, l_int_file);
    in.close();
}

TEST(SDSL, RemoveRange) {
    annotate::bv_t bv(128);
    for (size_t _i = 0; _i < -1llu - -1llu / 10; _i += -1llu / 10) {
        for (size_t _j = 0; _j < -1llu - -1llu / 10; _j += -1llu / 10) {
            bv.set_int(0, _i);
            bv.set_int(64, _j);
            for (size_t i = 0; i < 128; ++i) {
                for (size_t j = 1; j + i <= 128; ++j) {
                    annotate::bv_t bvv(128 - j);
                    if (!i) {
                        if (128 <= j + 64) {
                            bvv.set_int(0, bv.get_int(j, 128 - j), 128 - j);
                        } else {
                            bvv.set_int(0, bv.get_int(j));
                            bvv.set_int(64, bv.get_int(j + 64, 64 - j), 64 - j);
                        }
                    } else {
                        if (i >= 64) {
                            bvv.set_int(0, bv.get_int(0));
                            bvv.set_int(64, bv.get_int(64, i - 64), i - 64);
                            if (128 > i + j)
                                bvv.set_int(i, bv.get_int(i + j, 128 - i - j), 128 - i - j);
                        } else {
                            bvv.set_int(0, bv.get_int(0, i), i);
                            if (i + j <= 64) {
                                bvv.set_int(i, bv.get_int(i + j));
                                bvv.set_int(i + 64, bv.get_int(i + j + 64, 64 - i - j), 64 - i - j);
                            } else {
                                bvv.set_int(i, bv.get_int(i + j, 128 - i - j), 128 - i - j);
                            }
                        }
                    }
                    annotate::bv_t bvc(128);
                    bvc.set_int(0, bv.get_int(0));
                    bvc.set_int(64, bv.get_int(64));
                    ASSERT_EQ(bv, bvc);
                    annotate::bv_t bvd = annotate::remove_range(annotate::beta_t(bvc), i, i + j);
                    ASSERT_EQ(bvv, bvd) << i << "\t" << j;
                }
            }
        }
    }
}

TEST(SDSL, InsertZeros) {
    constexpr size_t size = 1000;
    for (size_t i = 0; i <= size; ++i) {
        annotate::bv_t bv(size, 1);
        annotate::bv_t bv_s(size + 1, 1);
        bv_s[i] = 0;
        ASSERT_EQ(bv_s, annotate::insert_zeros(bv, 1, i)) << i;
    }
}

TEST(SDSL, InsertRange) {
    constexpr size_t size = 1000;
    for (size_t i = 0; i < size + 2; ++i) { //insert size
        for (size_t j = 0; j <= size; ++j) { //insert position
            annotate::bv_t bv(size, 1);
            annotate::bv_t bv_s(size + i, 1);
            for (size_t k = j; k < j + i; ++k) {
                bv_s[k] = 0;
            }
            annotate::bv_t zeros(i);
            ASSERT_EQ(bv_s, annotate::insert_range(bv, zeros, j)) << j << " " << i;
        }
    }
}
