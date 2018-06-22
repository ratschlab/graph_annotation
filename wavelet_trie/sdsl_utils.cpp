#include "sdsl_utils.hpp"


namespace annotate {

constexpr size_t extract_bs_ = 64;

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
    for (; j + extract_bs_ <= i; j += extract_bs_) {
        merged.set_int(j, target.get_int(j, extract_bs_), extract_bs_);
    }
    if (i > j)
        merged.set_int(j, target.get_int(j, i - j), i - j);
    j = 0;
    for (; j + extract_bs_ <= source.size(); j += extract_bs_) {
        merged.set_int(i + j, source.get_int(j, extract_bs_), extract_bs_);
    }
    if (source.size() > j)
        merged.set_int(i + j, source.get_int(j, source.size() - j), source.size() - j);
    j = i;
    for (; j + extract_bs_ <= target.size(); j += extract_bs_) {
        merged.set_int(source.size() + j, target.get_int(j, extract_bs_), extract_bs_);
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

template <typename Vector>
bv_t remove_bits(const Vector &source, const std::vector<pos_t> &js) {
    if (js.size() > source.size()) {
        std::cerr << "too many indices\n";
        exit(1);
    }
    bv_t spliced;
    spliced.resize(source.size() - js.size());
    if (!spliced.size())
        return spliced;
    size_t j = 0;
    size_t i = 0;
    for (; i < js.size(); ++i) {
        assert(js[i] < source.size());
        for (; j + 64 <= js[i]; j += 64) {
            spliced.set_int(j - i, source.get_int(j));
        }
        if (js[i] > j) {
            spliced.set_int(j - i, source.get_int(j, js[i] - j), js[i] - j);
        }
        j = js[i] + 1;
    }
    for (; j + 64 <= source.size(); j += 64) {
        spliced.set_int(j - i, source.get_int(j));
    }
    if (j != source.size()) {
        spliced.set_int(j - i, source.get_int(j, source.size() - j), source.size() - j);
    }
    return spliced;
}
template bv_t remove_bits(const bv_t&, const std::vector<pos_t>&);
template bv_t remove_bits(const rrr_t&, const std::vector<pos_t>&);

template <typename Vector>
bv_t swap_bits(const Vector &source, const std::vector<std::pair<pos_t, pos_t>> &from_to) {
    bv_t swapped;
    swapped.resize(source.size());
    size_t j = 0;
    for (; j + 64 <= source.size(); j += 64) {
        swapped.set_int(j, source.get_int(j));
    }
    if (j != source.size()) {
        swapped.set_int(j, source.get_int(j, source.size() - j), source.size() - j);
    }
    for (auto &coords : from_to) {
        assert(coords.first < source.size());
        assert(coords.second < source.size());
        if (coords.first == coords.second)
            continue;
        bool temp = swapped[coords.second];
        swapped[coords.second] = swapped[coords.first];
        swapped[coords.first] = temp;
    }
    return swapped;
}
template bv_t swap_bits(const bv_t&, const std::vector<std::pair<pos_t, pos_t>>&);
template bv_t swap_bits(const rrr_t&, const std::vector<std::pair<pos_t, pos_t>>&);

template <class BitContainer>
void rearrange_bits(BitContainer &moved_c, const std::vector<std::pair<pos_t, pos_t>> &from_to) {
    size_t size = moved_c.size();
    auto it = from_to.begin();
    size_t to_offset = 0;
    size_t from_offset = 0;
    //size_t j = 0;

    //while (j < size && it != from_to.end()) {
    while (it != from_to.end()) {
        /*
        if (j + from_offset != it->first) {
            ++j;
            continue;
        }
        */
        //size_t j = it->first - from_offset;
        if (it->first < it->second) {
            auto jt = it;
            ++jt;
            while (jt != from_to.end()
                    && jt->second == it->second
                    && jt->first < jt->second
                    && jt->first == (jt - 1)->first + 1) {
                ++jt;
            }
            auto kt = moved_c.begin() + it->first - from_offset;
            std::rotate(kt,
                        kt + (jt - it),
                        moved_c.begin() + std::min(size, it->second + to_offset));
            //std::rotate(moved_c.begin() + j,
            //            moved_c.begin() + j + (jt - it),
            //            moved_c.begin() + std::min(size, it->second + to_offset));
            //if (it->second > j)
            if (it->second + from_offset > it->first)
                from_offset += jt - it;
                //++from_offset;
            it = jt;
            //continue;
        } else if (it->first > it->second) {
            auto jt = it;
            ++jt;
            while (jt != from_to.end()
                    && jt->second == it->second
                    && jt->first > jt->second
                    && jt->first == (jt - 1)->first + 1) {
                ++jt;
            }
            auto kt = moved_c.rend() - it->first + from_offset;
            std::rotate(kt - (jt - it),
                        kt,
                        moved_c.rend() - std::min(size, it->second + to_offset));
            //std::rotate(moved_c.rend() - j - (jt - it),
            //            moved_c.rend() - j,
            //            moved_c.rend() - std::min(size, it->second + to_offset));
            to_offset += jt - it;
            it = jt;
            //continue;
            //++to_offset;
        } else {
            ++it;
        }
        //++j;
        //}
        //++it;
    }
}

template <typename Vector>
std::vector<char> extract_bits(const Vector &source) {
    return std::vector<char>(source.begin(), source.end());
}

template <typename Vector>
bv_t copy_bits(const Vector &source) {
    bv_t moved;
    moved.resize(source.size());
    auto it = moved.begin();
    for (const auto &c : source) {
        *(it++) = c;
    }
    return moved;
}

template <typename Vector>
bv_t move_bits(const Vector &source,
               const std::vector<std::pair<pos_t, pos_t>> &from_to,
               pos_t remove_after) {
    //assert(std::is_sorted(from_to.begin(), from_to.end()));
    if (from_to.empty()) {
        return source;
    }

    //std::vector<char> moved_c(source.begin(), source.end());
    auto moved_c = extract_bits(source);
    rearrange_bits(moved_c, from_to);

    std::ignore = remove_after;
    auto moved = copy_bits(moved_c);
    /*
    bv_t moved;
    moved.resize(source.size());
    //moved.resize(std::min(remove_after, source.size()));
    //moved_c.resize(moved.size());
    auto _it = moved.begin();
    for (auto &c : moved_c) {
        *(_it++) = c;
    }
    */

    return moved;
}
template bv_t move_bits(const bv_t&, const std::vector<std::pair<pos_t, pos_t>>&, pos_t);
//template bv_t move_bits(const rrr_t&, const std::vector<std::pair<size_t, size_t>>&, size_t);


}; // annotate
