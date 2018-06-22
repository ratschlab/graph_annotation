#ifndef __SDSL_UTILS___
#define __SDSL_UTILS___

#include <sdsl/wavelet_trees.hpp>
#include <fstream>


namespace annotate {

typedef uint32_t pos_t;

constexpr size_t block_size_ = 255;
constexpr size_t sample_rate_ = 16;

typedef sdsl::bit_vector bv_t;
typedef sdsl::rrr_vector<block_size_, sdsl::int_vector<>, sample_rate_> rrr_t;

template <typename Vector>
bv_t insert_zeros(const Vector &target, const size_t count = 0, const size_t i = 0);

template <typename Vector1, typename Vector2>
bv_t insert_range(const Vector1 &target, const Vector2 &source, const size_t i = 0);

template <typename Vector>
bv_t remove_range(const Vector &source, const size_t begin, const size_t end);
template <typename Vector>
bv_t remove_bits(const Vector &source, const std::vector<pos_t> &js);

template <typename Vector>
bv_t swap_bits(const Vector &source, const std::vector<std::pair<pos_t, pos_t>> &from_to);

template <typename Vector>
bv_t move_bits(const Vector &source,
               const std::vector<std::pair<pos_t, pos_t>> &from_to,
                pos_t remove_after = static_cast<pos_t>(-1));

}; // annotate

#endif // __SDSL_UTILS___
