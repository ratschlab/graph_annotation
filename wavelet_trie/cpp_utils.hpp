#ifndef __CPP_UTILS___
#define __CPP_UTILS___

#include <iostream>
#include <set>
#include <vector>
#include <boost/multiprecision/gmp.hpp>


namespace annotate {

typedef boost::multiprecision::mpz_int cpp_int;

bool is_nonzero(const cpp_int &a);
bool is_nonzero(const std::set<size_t> &a);
bool is_nonzero(const std::vector<size_t> &a);

bool bit_test(const cpp_int &a, const size_t &col);
bool bit_test(const std::set<size_t> &a, const size_t &col);
bool bit_test(const std::vector<size_t> &a, const size_t &col);

void bit_set(mpz_t &a_d, const size_t &col);
void bit_set(cpp_int &a, const size_t &col);
void bit_set(std::set<size_t> &a, const size_t &col);
void bit_set(std::vector<size_t> &a, const size_t &col);

void bit_unset(cpp_int &a, const size_t &col);
void bit_unset(std::set<size_t> &a, const size_t &col);
void bit_unset(std::vector<size_t> &a, const size_t &col);

template <class Container>
void bit_toggle(Container &a, const size_t &col);

size_t next_bit(const cpp_int &a, const size_t &col);
size_t next_bit(const std::set<size_t> &a, const size_t &col);
size_t next_bit(const std::vector<size_t> &a, const size_t &col);

void clear_after(mpz_t &a_d, const size_t &col);
void clear_after(cpp_int &a, const size_t &col);
void clear_after(std::set<size_t> &a, const size_t &col);
void clear_after(std::vector<size_t> &a, const size_t &col);

size_t msb(const cpp_int &a);
size_t msb(const std::set<size_t> &a);
size_t msb(const std::vector<size_t> &a);

size_t lsb(const cpp_int &a);
size_t lsb(const std::set<size_t> &a);
size_t lsb(const std::vector<size_t> &a);

size_t popcount(const cpp_int &a);

size_t serialize(std::ostream &out, const cpp_int &l_int);

cpp_int load(std::istream &in);

cpp_int pack_indices(const std::vector<size_t> &indices);

}; // annotate

#endif // __CPP_UTILS___
