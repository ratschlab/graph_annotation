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

std::vector<std::vector<size_t>> generate_indices(std::vector<std::set<size_t>> &bits) {
    std::vector<std::vector<size_t>> nums(bits.size());
    for (size_t i = 0; i < nums.size(); ++i) {
        nums[i].insert(nums[i].end(), bits[i].begin(), bits[i].end());
    }
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
template <class Container>
void generate_wtr(std::vector<annotate::WaveletTrie> &wtrs, Container nums, size_t step = -1llu) {
    auto it = nums.begin();
    while (step < -1llu && it + step <= nums.end()) {
        wtrs.emplace_back(it, it + step);
        it += step;
    }
    if (it != nums.end())
        wtrs.emplace_back(it, nums.end());
}

annotate::WaveletTrie merge_wtrs(std::vector<annotate::WaveletTrie> &wtrs) {
    annotate::WaveletTrie wts;
    for (auto &wtr : wtrs) {
        wts.insert(std::move(wtr));
    }
    return wts;
}

annotate::WaveletTrie merge_wtrs_copy(std::vector<annotate::WaveletTrie> &wtrs) {
    annotate::WaveletTrie wts;
    for (const auto &wtr : wtrs) {
        wts.insert(wtr);
    }
    return wts;
}

annotate::WaveletTrie test_wtr(size_t i, size_t p, int mode = 0) {
    std::vector<annotate::WaveletTrie> wtrs;
    auto nums = generate_nums(bits[i]);
    if (mode == 0) {
        generate_wtr(wtrs, nums);
    } else if (mode == 1) {
        auto nums = generate_indices(bits[i]);
        generate_wtr(wtrs, nums);
    } else {
        generate_wtr(wtrs, bits[i]);
    }
    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, nums, std::to_string(i));
    return wts;
}
annotate::WaveletTrie test_wtr_step(size_t i, size_t p, int mode = 0) {
    std::vector<annotate::WaveletTrie> wtrs;
    auto nums = generate_nums(bits[i]);
    constexpr size_t step = 2;

    if (mode == 0) {
        generate_wtr(wtrs, nums, step);
    } else if (mode == 1) {
        auto nums = generate_indices(bits[i]);
        generate_wtr(wtrs, nums, step);
    } else {
        generate_wtr(wtrs, bits[i], step);
    }

    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, nums);
    return wts;
}

annotate::WaveletTrie test_wtr_copy(size_t i, size_t p, int mode = 0) {
    std::vector<annotate::WaveletTrie> wtrs;
    auto nums = generate_nums(bits[i]);

    if (mode == 0) {
        generate_wtr(wtrs, nums);
    } else if (mode == 1) {
        auto nums = generate_indices(bits[i]);
        generate_wtr(wtrs, nums);
    } else {
        generate_wtr(wtrs, bits[i]);
    }

    auto wts = merge_wtrs_copy(wtrs);
    check_wtr(wts, nums, std::to_string(i));
    return wts;
}
annotate::WaveletTrie test_wtr_copy_step(size_t i, size_t p, int mode = 0) {
    std::vector<annotate::WaveletTrie> wtrs;
    auto nums = generate_nums(bits[i]);
    constexpr size_t step = 2;

    if (mode == 0) {
        generate_wtr(wtrs, nums, step);
    } else if (mode == 1) {
        auto nums = generate_indices(bits[i]);
        generate_wtr(wtrs, nums, step);
    } else {
        generate_wtr(wtrs, bits[i], step);
    }
    auto wts = merge_wtrs_copy(wtrs);
    check_wtr(wts, nums);
    return wts;
}

annotate::WaveletTrie test_wtr_pairs(size_t i, size_t j, size_t p, int mode = 0) {
    auto nums1 = generate_nums(bits[i]);
    auto nums2 = generate_nums(bits[j]);
    std::vector<annotate::WaveletTrie> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    if (mode == 0) {
        generate_wtr(wtrs, nums1);
        generate_wtr(wtrs, nums2);
    } else if (mode == 1) {
        auto nums1 = generate_indices(bits[i]);
        auto nums2 = generate_indices(bits[j]);
        generate_wtr(wtrs, nums1);
        generate_wtr(wtrs, nums2);
    } else if (mode == 2) {
        generate_wtr(wtrs, bits[i]);
        generate_wtr(wtrs, bits[j]);
    }
    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j));
    return wts;
}
annotate::WaveletTrie test_wtr_pairs_step(size_t i, size_t j, size_t p, int mode = 0) {
    auto nums1 = generate_nums(bits[i]);
    auto nums2 = generate_nums(bits[j]);
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    if (mode == 0) {
        generate_wtr(wtrs, nums1, step);
        generate_wtr(wtrs, nums2, step);
    } else if (mode == 1) {
        auto nums1 = generate_indices(bits[i]);
        auto nums2 = generate_indices(bits[j]);
        generate_wtr(wtrs, nums1, step);
        generate_wtr(wtrs, nums2, step);
    } else {
        generate_wtr(wtrs, bits[i], step);
        generate_wtr(wtrs, bits[j], step);
    }
    auto wts = merge_wtrs(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(step));
    return wts;
}

annotate::WaveletTrie test_wtr_pairs_copy(size_t i, size_t j, size_t p, int mode = 0) {
    auto nums1 = generate_nums(bits[i]);
    auto nums2 = generate_nums(bits[j]);
    std::vector<annotate::WaveletTrie> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    if (mode == 0) {
        generate_wtr(wtrs, nums1);
        generate_wtr(wtrs, nums2);
    } else if (mode == 1) {
        auto nums1 = generate_indices(bits[i]);
        auto nums2 = generate_indices(bits[j]);
        generate_wtr(wtrs, nums1);
        generate_wtr(wtrs, nums2);
    } else {
        generate_wtr(wtrs, bits[i]);
        generate_wtr(wtrs, bits[j]);
    }
    auto wts = merge_wtrs_copy(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j));
    return wts;
}

annotate::WaveletTrie test_wtr_pairs_copy_step(size_t i, size_t j, size_t p, int mode = 0) {
    auto nums1 = generate_nums(bits[i]);
    auto nums2 = generate_nums(bits[j]);
    constexpr size_t step = 2;
    std::vector<annotate::WaveletTrie> wtrs;
    std::vector<annotate::cpp_int> ref;
    ref.reserve(nums1.size() + nums2.size());
    ref.insert(ref.end(), nums1.begin(), nums1.end());
    ref.insert(ref.end(), nums2.begin(), nums2.end());

    if (mode == 0) {
        generate_wtr(wtrs, nums1, step);
        generate_wtr(wtrs, nums2, step);
    } else if (mode == 1) {
        auto nums1 = generate_indices(bits[i]);
        auto nums2 = generate_indices(bits[j]);
        generate_wtr(wtrs, nums1, step);
        generate_wtr(wtrs, nums2, step);
    } else { 
        generate_wtr(wtrs, bits[i], step);
        generate_wtr(wtrs, bits[j], step);
    }
    auto wts = merge_wtrs_copy(wtrs);
    check_wtr(wts, ref, std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(step));
    return wts;
}

void check_wtr_vector(std::vector<annotate::WaveletTrie> &wtrs) {
    for (size_t i = 1; i < wtrs.size(); ++i) {
        ASSERT_EQ(wtrs[i - 1], wtrs[i]) << "EQUAL FAIL: " << i << std::endl;
    }
}

std::vector<size_t> num_threads = {1, 4};

TEST(WaveletTrie, TestSingle) {
    std::vector<annotate::WaveletTrie> wtrs;
    for (size_t i = 0; i < bits.size(); ++i) {
        wtrs.clear();
        for (auto p : num_threads) {
            wtrs.push_back(test_wtr(i, p, 0));
            wtrs.push_back(test_wtr(i, p, 1));
            wtrs.push_back(test_wtr(i, p, 2));
            wtrs.push_back(test_wtr_step(i, p, 0));
            wtrs.push_back(test_wtr_step(i, p, 1));
            wtrs.push_back(test_wtr_step(i, p, 2));
            wtrs.push_back(test_wtr_copy(i, p, 0));
            wtrs.push_back(test_wtr_copy(i, p, 1));
            wtrs.push_back(test_wtr_copy(i, p, 2));
            wtrs.push_back(test_wtr_copy_step(i, p, 0));
            wtrs.push_back(test_wtr_copy_step(i, p, 1));
            wtrs.push_back(test_wtr_copy_step(i, p, 2));
        }
        check_wtr_vector(wtrs);
    }
}

TEST(WaveletTrie, TestPairs) {
    std::vector<annotate::WaveletTrie> wtrs;
    for (size_t i = 0; i < bits.size(); ++i) {
        for (size_t j = 0; j < bits.size(); ++j) {
            wtrs.clear();
            for (auto p : num_threads) {
                wtrs.push_back(test_wtr_pairs(i, j, p, 0));
                wtrs.push_back(test_wtr_pairs(i, j, p, 1));
                wtrs.push_back(test_wtr_pairs(i, j, p, 2));
                wtrs.push_back(test_wtr_pairs_step(i, j, p, 0));
                wtrs.push_back(test_wtr_pairs_step(i, j, p, 1));
                wtrs.push_back(test_wtr_pairs_step(i, j, p, 2));
                wtrs.push_back(test_wtr_pairs_copy(i, j, p, 0));
                wtrs.push_back(test_wtr_pairs_copy(i, j, p, 1));
                wtrs.push_back(test_wtr_pairs_copy(i, j, p, 2));
                wtrs.push_back(test_wtr_pairs_copy_step(i, j, p, 0));
                wtrs.push_back(test_wtr_pairs_copy_step(i, j, p, 1));
                wtrs.push_back(test_wtr_pairs_copy_step(i, j, p, 2));
            }
            check_wtr_vector(wtrs);
        }
    }
}


