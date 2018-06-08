#ifndef __SDSL_UTILS___
#define __SDSL_UTILS___

#include <sdsl/wavelet_trees.hpp>


namespace annotate {

typedef sdsl::bit_vector bv_t;
typedef sdsl::rrr_vector<255,sdsl::int_vector<>,16> rrr_t;

template <typename Vector>
bv_t insert_zeros(const Vector &target, const size_t count = 0, const size_t i = 0);

template <typename Vector1, typename Vector2>
bv_t insert_range(const Vector1 &target, const Vector2 &source, const size_t i = 0);

template <typename Vector>
bv_t remove_range(const Vector &source, const size_t begin, const size_t end);

}; // annotate

#endif // __SDSL_UTILS___
