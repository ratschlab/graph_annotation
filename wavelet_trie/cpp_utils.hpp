#ifndef __CPP_UTILS___
#define __CPP_UTILS___

#include <iostream>
#include <set>
#include <vector>
#include <boost/multiprecision/gmp.hpp>


namespace annotate {

typedef boost::multiprecision::mpz_int cpp_int;
typedef uint32_t pos_t;

bool is_nonzero(const cpp_int &a);
bool is_nonzero(const std::set<pos_t> &a);
bool is_nonzero(const std::vector<pos_t> &a);

bool bit_test(const cpp_int &a, size_t col);
bool bit_test(const std::set<pos_t> &a, size_t col);
bool bit_test(const std::vector<pos_t> &a, size_t col);

void bit_set(mpz_t &a_d, size_t col);
void bit_set(cpp_int &a, size_t col);
void bit_set(std::set<pos_t> &a, size_t col);
void bit_set(std::vector<pos_t> &a, size_t col);

void bit_unset(cpp_int &a, size_t col);
void bit_unset(std::set<pos_t> &a, size_t col);
void bit_unset(std::vector<pos_t> &a, size_t col);

template <class Container>
void bit_toggle(Container &a, size_t col);

pos_t next_bit(const cpp_int &a, size_t col);
pos_t next_bit(const std::set<pos_t> &a, size_t col);
pos_t next_bit(const std::vector<pos_t> &a, size_t col);

void clear_after(mpz_t &a_d, size_t col);
void clear_after(cpp_int &a, size_t col);
void clear_after(std::set<pos_t> &a, size_t col);
void clear_after(std::vector<pos_t> &a, size_t col);

pos_t msb(const cpp_int &a);
pos_t msb(const std::set<pos_t> &a);
pos_t msb(const std::vector<pos_t> &a);

pos_t lsb(const cpp_int &a);
pos_t lsb(const std::set<pos_t> &a);
pos_t lsb(const std::vector<pos_t> &a);

size_t popcount(const cpp_int &a);

size_t serialize(std::ostream &out, const cpp_int &l_int);

cpp_int load(std::istream &in);

template <class Container>
cpp_int pack_indices(const Container &indices);

}; // annotate

#endif // __CPP_UTILS___
