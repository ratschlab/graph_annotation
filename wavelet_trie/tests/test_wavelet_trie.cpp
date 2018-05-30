#include "gtest/gtest.h"
#include "wavelet_trie.hpp"

std::vector<annotate::cpp_int> generate_nums(std::vector<std::set<size_t>> &bits) {
    std::vector<annotate::cpp_int> nums;
    for (auto &num : bits) {
        nums.emplace_back(0);
        for (auto &i : num) {
            annotate::bit_set(nums.back(), i);
        }
    }
    return nums;
}

void test_wtr(std::vector<annotate::cpp_int>&& nums) {
    std::vector<annotate::cpp_int> ref(nums.begin(), nums.end());
    annotate::WaveletTrie wt(nums);
    ASSERT_EQ(ref.size(), wt.size());
    for (size_t i = 0; i < ref.size(); ++i) {
        EXPECT_EQ(ref.at(i), wt.at(i));
    }
}

TEST(WaveletTrie, Test0) {
    std::vector<std::set<size_t>> bits;
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test1) {
    std::vector<std::set<size_t>> bits {
        {4},
        {2},
        {4},
        {2},
        {4},
        {2},
        {2}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test2) {
    std::vector<std::set<size_t>> bits {
        {3},
        {2},
        {2, 3}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test3) {
    std::vector<std::set<size_t>> bits {
        {3},
        {2},
        {2, 3},
        {3, 4, 5}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test4) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 3}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test5) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 2}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test6) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 3, 4}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test7) {
    std::vector<std::set<size_t>> bits {
        {1, 4},
        {1, 3, 4}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test8) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 7},
        {1, 5, 7},
        {1, 3, 7}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test9) {
    std::vector<std::set<size_t>> bits {
        {0},
        {1},
        {1}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test10) {
    std::vector<std::set<size_t>> bits {
        {1},
        {5}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test11) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 5}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test12) {
    std::vector<std::set<size_t>> bits {
        {0},
        {1},
        {5}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test13) {
    std::vector<std::set<size_t>> bits {
        {0},
        {1},
        {5},
        {1},
        {5},
        {1},
        {5}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test14) {
    std::vector<std::set<size_t>> bits {
        {1},
        {3},
        {1, 3},
        {1}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test15) {
    std::vector<std::set<size_t>> bits {
        {1, 3}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test16) {
    std::vector<std::set<size_t>> bits {
        {1},
        {1, 3},
        {1, 3, 4},
        {1},
        {1, 3},
        {1, 3}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test17) {
    std::vector<std::set<size_t>> bits {
        {3, 4},
        {3, 4, 6},
        {3, 4, 6, 7, 8, 9}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test18) {
    std::vector<std::set<size_t>> bits {
        {3, 4, 6, 7, 8, 9, 10},
        {3, 4, 6, 7, 9}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test19) {
    std::vector<std::set<size_t>> bits {
        {1},
        {3},
        {1},
        {3},
        {1},
        {3},
        {1, 3}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test20) {
    std::vector<std::set<size_t>> bits {
        {1, 3, 5},
        {1, 3, 4, 6},
        {1, 3, 4},
        {1, 2, 3, 5},
        {1, 2, 3, 4, 6},
        {1, 2, 3, 4}
    };
    test_wtr(generate_nums(bits));
}

TEST(WaveletTrie, Test21) {
    std::vector<std::set<size_t>> bits {
        {1, 3, 5},
        {1, 3, 4, 6},
        {1, 3, 4},
        {0, 1, 2, 3, 5},
        {0, 1, 2, 3, 4, 6},
        {1, 2, 3, 4, 4}
    };
    test_wtr(generate_nums(bits));
}

