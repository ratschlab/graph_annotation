#include "gtest/gtest.h"
#include "wavelet_trie.hpp"


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


std::vector<annotate::cpp_int> generate_nums(std::vector<std::set<size_t>> &bits) {
    std::vector<annotate::cpp_int> nums;
    nums.reserve(bits.size());
    std::transform(bits.begin(), bits.end(), std::back_inserter(nums), annotate::pack_indices);
    return nums;
}

void check_wtr(annotate::WaveletTrie &wt,
               const std::vector<annotate::cpp_int> &nums,
               const std::string &message = std::string()) {
    ASSERT_EQ(nums.size(), wt.size()) << message << ": size fail" << std::endl;
    for (size_t i = 0; i < nums.size(); ++i) {
        EXPECT_EQ(nums.at(i), wt.at(i)) << message << ":" << std::to_string(i) << std::endl;
    }
}

// Note: the vector needs to be a copy, since WaveletTrie modifies it
void generate_wtr(std::vector<annotate::WaveletTrie> &wtrs, std::vector<annotate::cpp_int> nums, size_t step = -1llu) {
    auto it = nums.begin();
    while (step < -1llu && it + step <= nums.end()) {
        wtrs.emplace_back(it, it + step);
        it += step;
    }
    if (it != nums.end())
        wtrs.emplace_back(it, nums.end());
}

// TODO: this only works as a pointer, since destructor for int_vector fails when static
annotate::WaveletTrie merge_wtrs(std::vector<annotate::WaveletTrie> &wtrs) {
    annotate::WaveletTrie wts;
    for (auto &wtr : wtrs) {
        wts.insert(std::move(wtr));
    }
    return wts;
}

void test_wtr(size_t i) {
    auto nums = generate_nums(bits[i]);
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie> wtrs;

    generate_wtr(wtrs, nums, step);
    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, nums);
    wtrs.clear();

    generate_wtr(wtrs, nums);
    wts = merge_wtrs(wtrs);
    check_wtr(wts, nums, std::to_string(i));
}

void test_wtr_pairs(size_t i, size_t j) {
    auto nums1 = generate_nums(bits[i]);
    auto nums2 = generate_nums(bits[j]);
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    generate_wtr(wtrs, nums1, step);
    generate_wtr(wtrs, nums2, step);
    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(step));
    wtrs.clear();

    generate_wtr(wtrs, nums1);
    generate_wtr(wtrs, nums2);
    wts = merge_wtrs(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j));
    wtrs.clear();
}

TEST(WaveletTrie, TestSingle) {
    for (size_t i = 0; i < bits.size(); ++i) {
        test_wtr(i);
    }
}

TEST(WaveletTrie, TestPairs) {
    for (size_t i = 0; i < bits.size(); ++i) {
        for (size_t j = 0; j < bits.size(); ++j) {
            test_wtr_pairs(i, j);
        }
    }
}

