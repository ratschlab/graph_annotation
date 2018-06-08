#include "sdsl_utils.hpp"


namespace annotate {

//TODO: align get_int to blocks in source
template <typename Vector>
bv_t insert_zeros(const Vector &target, const size_t count, const size_t i) {
    //if (!count) {
    //    return target;
    //}
    if (!target.size()) {
        return bv_t(count);
    }
    bv_t merged;
    size_t j;
    //assume limb size of 64
    if (i) {
        //merged = target;
        merged.resize(target.size() + count);
        j = 0;
        for (; j + 64 <= i; j += 64) {
            merged.set_int(j, target.get_int(j));
        }
        if (i > j)
            merged.set_int(j, target.get_int(j, i - j), i - j);

        if (count) {
            j = ((i + 63) & -64llu);
            //WARNING: the region from merged.size() to j is not initialized
            merged.set_int(i, 0, std::min(j, static_cast<size_t>(merged.size())) - i);
            std::fill(merged.data() + (j >> 6), merged.data() + (merged.capacity() >> 6), 0);
        }
    } else {
        merged = bv_t(count);
        merged.resize(target.size() + count);
    }
    j = i;
    for (; j + 64 <= target.size(); j += 64) {
        merged.set_int(count + j, target.get_int(j));
    }
    if (target.size() > j)
        merged.set_int(count + j, target.get_int(j, target.size() - j), target.size() - j);
    return merged;
}
template bv_t insert_zeros(const bv_t&,  const size_t, const size_t);
template bv_t insert_zeros(const rrr_t&, const size_t, const size_t);

template <typename Vector1, typename Vector2>
bv_t insert_range(const Vector1 &target, const Vector2 &source, const size_t i) {
    bv_t merged;
    merged.resize(target.size() + source.size());
    size_t j = 0;
    for (; j + 64 <= i; j += 64) {
        merged.set_int(j, target.get_int(j));
    }
    if (i > j)
        merged.set_int(j, target.get_int(j, i - j), i - j);
    j = 0;
    for (; j + 64 <= source.size(); j += 64) {
        merged.set_int(i + j, source.get_int(j));
    }
    if (source.size() > j)
        merged.set_int(i + j, source.get_int(j, source.size() - j), source.size() - j);
    j = i;
    for (; j + 64 <= target.size(); j += 64) {
        merged.set_int(source.size() + j, target.get_int(j));
    }
    if (target.size() > j)
        merged.set_int(source.size() + j, target.get_int(j, target.size() - j), target.size() - j);
    return merged;
}
template bv_t insert_range(const bv_t&,  const bv_t&,  const size_t);
template bv_t insert_range(const rrr_t&, const rrr_t&, const size_t);
template bv_t insert_range(const bv_t&,  const rrr_t&, const size_t);
template bv_t insert_range(const rrr_t&, const bv_t&,  const size_t);

template <typename Vector>
bv_t remove_range(const Vector &source, const size_t begin, const size_t end) {
    if (begin > end) {
        std::cerr << "begin > end\n";
        exit(1);
    }
    if (begin >= source.size()) {
        std::cerr << "begin >= source.size()\n";
        exit(1);
    }
    if (end > source.size()) {
        std::cerr << "end > source.size()\n";
        exit(1);
    }

    bv_t spliced;
    spliced.resize(source.size() + begin - end);
    size_t j = 0;

    //first chunk
    for (; j + 64 <= begin; j += 64) {
        spliced.set_int(j, source.get_int(j));
    }
    if (begin > j) {
        spliced.set_int(j, source.get_int(j, begin - j), begin - j);
        j = begin;
    }
    //TODO: align this to limbs
    size_t i = end;
    for (; i + 64 <= source.size(); i += 64) {
        spliced.set_int(j, source.get_int(i));
        j += 64;
    }
    if (i != source.size()) {
        spliced.set_int(j, source.get_int(i, source.size() - i), source.size() - i);
    }
    return spliced;
}
template bv_t remove_range(const bv_t&,  const size_t, const size_t);
template bv_t remove_range(const rrr_t&, const size_t, const size_t);

}; // annotate
