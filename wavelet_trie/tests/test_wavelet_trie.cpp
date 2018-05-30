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

void check_wtr(annotate::WaveletTrie &wt, const std::vector<annotate::cpp_int> &nums) {
    ASSERT_EQ(nums.size(), wt.size());
    for (size_t i = 0; i < nums.size(); ++i) {
        EXPECT_EQ(nums.at(i), wt.at(i));
    }
}

// Note: the vector needs to be a copy, since WaveletTrie modifies it
void generate_wtr(std::vector<annotate::WaveletTrie*> &wtrs, std::vector<annotate::cpp_int> nums, size_t step = -1llu) {
    auto it = nums.begin();
    while (step < -1llu && it + step <= nums.end()) {
        wtrs.push_back(new annotate::WaveletTrie(it, it + step));
        it += step;
    }
    if (it != nums.end())
        wtrs.push_back(new annotate::WaveletTrie(it, nums.end()));
}

// TODO: this only works as a pointer, since destructor for int_vector fails when static
annotate::WaveletTrie* merge_wtrs(std::vector<annotate::WaveletTrie*> &wtrs) {
    annotate::WaveletTrie *wts = new annotate::WaveletTrie();
    for (auto &wtr : wtrs) {
        wts->insert(*wtr);
        delete wtr;
    }
    return wts;
}

void test_wtr(std::vector<annotate::cpp_int>&& nums) {
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie*> wtrs;

    generate_wtr(wtrs, nums, step);
    auto wts = merge_wtrs(wtrs);
    check_wtr(*wts, nums);
    delete wts;
    wtrs.clear();

    generate_wtr(wtrs, nums);
    wts = merge_wtrs(wtrs);
    check_wtr(*wts, nums);
    delete wts;
}

void test_wtr_pairs(std::vector<annotate::cpp_int>&& nums1, std::vector<annotate::cpp_int>&& nums2) {
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie*> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    generate_wtr(wtrs, nums1, step);
    generate_wtr(wtrs, nums2, step);
    auto wts = merge_wtrs(wtrs);
    check_wtr(*wts, ref);
    delete wts;
    wtrs.clear();
}

std::vector<std::vector<std::set<size_t>>> bits {
    {}, // 0
    {   // 1
        {4},
        {2},
        {4},
        {2},
        {4},
        {2},
        {2}
    },
    {   // 2
        {3},
        {2},
        {2, 3}
    },
    {   // 3
        {3},
        {2},
        {2, 3},
        {3, 4, 5}
    },
    {
        {1},
        {1, 3}
    },
    {
        {1},
        {1, 2}
    },
    {
        {1},
        {1, 3, 4}
    },
    {
        {1, 4},
        {1, 3, 4}
    },
    {
        {1},
        {1, 7},
        {1, 5, 7},
        {1, 3, 7}
    },
    {
        {0},
        {1},
        {1}
    },
    {
        {1},
        {5}
    },
    {
        {1},
        {1, 5}
    },
    {
        {0},
        {1},
        {5}
    },
    {
        {0},
        {1},
        {5},
        {1},
        {5},
        {1},
        {5}
    },
    {
        {1},
        {3},
        {1, 3},
        {1}
    },
    {
        {1, 3}
    },
    {
        {1},
        {1, 3},
        {1, 3, 4},
        {1},
        {1, 3},
        {1, 3}
    },
    {
        {3, 4},
        {3, 4, 6},
        {3, 4, 6, 7, 8, 9}
    },
    {
        {3, 4, 6, 7, 8, 9, 10},
        {3, 4, 6, 7, 9}
    },
    {
        {1},
        {3},
        {1},
        {3},
        {1},
        {3},
        {1, 3}
    },
    {
        {1, 3, 5},
        {1, 3, 4, 6},
        {1, 3, 4},
        {1, 2, 3, 5},
        {1, 2, 3, 4, 6},
        {1, 2, 3, 4}
    },
    {
        {1, 3, 5},
        {1, 3, 4, 6},
        {1, 3, 4},
        {0, 1, 2, 3, 5},
        {0, 1, 2, 3, 4, 6},
        {1, 2, 3, 4, 4}
    }
};

TEST(WaveletTrie, TestSingle) {
    for (auto &bitset : bits) {
        test_wtr(generate_nums(bitset));
    }
}

TEST(WaveletTrie, TestPairs) {
    for (auto &bitset1 : bits) {
        for (auto &bitset2 : bits) {
            test_wtr_pairs(generate_nums(bitset1), generate_nums(bitset2));
        }
    }
}

