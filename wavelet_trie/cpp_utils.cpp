#include "cpp_utils.hpp"

namespace annotate {

bool is_nonzero(const cpp_int &a) {
    return a != 0;
}

bool is_nonzero(const std::set<size_t> &a) {
    return a.size();
}

bool is_nonzero(const std::vector<size_t> &a) {
    return a.size();
}

bool bit_test(const cpp_int &a, size_t col) {
    assert(col < -1llu);
    return mpz_tstbit(a.backend().data(), col);
}

bool bit_test(const std::set<size_t> &a, size_t col) {
    return a.find(col) != a.end();
}

bool bit_test(const std::vector<size_t> &a, size_t col) {
    return std::find(a.begin(), a.end(), col) != a.end();
}

void bit_set(mpz_t &a_d, size_t col) {
    assert(col < -1llu);
    mpz_setbit(a_d, col);
}

void bit_set(cpp_int &a, size_t col) {
    assert(col < -1llu);
    mpz_setbit(a.backend().data(), col);
}

void bit_set(std::set<size_t> &a, size_t col) {
    assert(col < -1llu);
    a.insert(col);
}

void bit_set(std::vector<size_t> &a, size_t col) {
    assert(col < -1llu);
    if (!bit_test(a, col)) {
        a.push_back(col);
    }
}
void bit_unset(cpp_int &a, size_t col) {
    assert(col < -1llu);
    mpz_clrbit(a.backend().data(), col);
}

void bit_unset(std::set<size_t> &a, size_t col) {
    auto index = a.find(col);
    if (index != a.end())
        a.erase(index);
}

void bit_unset(std::vector<size_t> &a, size_t col) {
    auto index = std::find(a.begin(), a.end(), col);
    if (index != a.end())
        a.erase(index);
}

template <class Container>
void bit_toggle(Container &a, size_t col) {
    assert(col < -1llu);
    if (bit_test(a, col))
        bit_unset(a, col);
    else
        bit_set(a, col);
}
template void bit_toggle(cpp_int&, size_t);
template void bit_toggle(std::vector<size_t>&, size_t);
template void bit_toggle(std::set<size_t>&, size_t);

size_t next_bit(const cpp_int &a, size_t col) {
    assert(col < -1llu);
    return mpz_scan1(a.backend().data(), col);
}

size_t next_bit(const std::set<size_t> &a, size_t col) {
    assert(col < -1llu);
    auto lbound = a.lower_bound(col);
    if (lbound == a.end())
        return -1llu;
    return *lbound;
}

size_t next_bit(const std::vector<size_t> &a, size_t col) {
    assert(col < -1llu);
    size_t lbound = -1llu;
    for (auto &index : a) {
        if (index >= col) {
            lbound = std::min(lbound, index);
            if (index == col)
                break;
        }
    }
    return lbound;
}

void clear_after(mpz_t &a_d, size_t col) {
    //assert(col < -1llu);
    if (col == -1llu)
        return;
    if (!col) {
        mpz_clear(a_d);
        mpz_init(a_d);
    } else {
        mpz_tdiv_r_2exp(a_d, a_d, col);
    }
}

void clear_after(cpp_int &a, size_t col) {
    assert(col < -1llu);
    mpz_t& a_d = a.backend().data();
    clear_after(a_d, col);
}

void clear_after(std::set<size_t> &a, size_t col) {
    a.erase(a.lower_bound(col), a.end());
}

void clear_after(std::vector<size_t> &a, size_t col) {
    auto it = a.begin();
    while (it != a.end()) {
        if (*it >= col) {
            a.erase(it);
        } else {
            ++it;
        }
    }
}

size_t msb(const cpp_int &a) {
    assert(a != 0);
    const mpz_t& a_mpz = a.backend().data();
    size_t i = mpz_scan1(a_mpz, 0);
    do {
        size_t j = mpz_scan1(a_mpz, i + 1);
        if (j == -1llu)
            break;
        i = j;
    } while (true);
    return i;
}

size_t msb(const std::set<size_t> &a) {
    assert(a.size());
    return *a.rbegin();
}

size_t msb(const std::vector<size_t> &a) {
    assert(a.size());
    return *std::max_element(a.begin(), a.end());
}

size_t lsb(const cpp_int &a) {
    assert(a != 0);
    return next_bit(a, 0);
}

size_t lsb(const std::set<size_t> &a) {
    assert(a.size());
    return *a.begin();
}

size_t lsb(const std::vector<size_t> &a) {
    assert(a.size());
    return *std::min_element(a.begin(), a.end());
}

size_t popcount(const cpp_int &a) {
    return mpz_popcount(a.backend().data());
}

size_t serialize(std::ostream &out, const cpp_int &l_int) {
    size_t a;
    void *l_int_raw = mpz_export(NULL, &a, 1, 1, 0, 0, l_int.backend().data());
    out.write(reinterpret_cast<char*>(&a), sizeof(a));
    out.write((char*)l_int_raw, a);
    free(l_int_raw);
    return a;
}

cpp_int load(std::istream &in) {
    cpp_int curint = 0;
    size_t a = 0;
    in.read(reinterpret_cast<char*>(&a), sizeof(a));
    void *l_int_raw = malloc(a);
    in.read((char*)l_int_raw, a);
    mpz_import(curint.backend().data(), a, 1, 1, 0, 0, l_int_raw);
    free(l_int_raw);
    return curint;
}

cpp_int pack_indices(const std::vector<size_t> &indices) {
    cpp_int num = 0;
    for (auto &i : indices) {
        bit_set(num, i);
    }
    return num;
}

}; // annotate
