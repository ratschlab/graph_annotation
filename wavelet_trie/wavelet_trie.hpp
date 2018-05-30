#ifndef __WAVELET_TRIE___
#define __WAVELET_TRIE___

#include <iostream>
#include <set>
#include <fstream>
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>
#include <vector>
#include <future>
#include <functional>
#include <boost/multiprecision/gmp.hpp>

#include "thread_pool.hpp"


namespace annotate {
    typedef boost::multiprecision::mpz_int cpp_int;
    typedef cpp_int annot_t;
    typedef sdsl::bit_vector bv_t;
    typedef sdsl::rrr_vector<255,sdsl::int_vector<>,16> rrr_t;
    typedef sdsl::sd_vector<> sd_t;

    typedef cpp_int alpha_t;
    typedef rrr_t beta_t;
    //typedef bv_t beta_t;
    typedef beta_t::rank_1_type rank1_t;

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

    size_t serialize(std::ostream &out, const cpp_int &l_int);

    cpp_int load(std::istream &in);

    cpp_int pack_indices(std::set<size_t> &indices);

    //template <typename Vector>
    //bv_t insert_zeros(const Vector &target, const size_t count = 0, const size_t i = 0);

    //template <typename Vector>
    //bv_t insert_range(const Vector &target, const Vector &source, const size_t i = 0);

    class Prefix {
        public:
            Prefix() { }
            Prefix(size_t col, bool allequal) : col(col), allequal(allequal) { }
            size_t col = -1llu;
            bool allequal = true;
    };

    class WaveletTrie {
        public:

            WaveletTrie();

            // copy constructor
            WaveletTrie(const WaveletTrie &other);

            // move constructor
            WaveletTrie(WaveletTrie&& other) noexcept;

            // assign
            WaveletTrie& operator=(const WaveletTrie& other);
            WaveletTrie& operator=(WaveletTrie&& other) noexcept;

            template <class Iterator>
            WaveletTrie(Iterator row_begin, Iterator row_end);
            //WaveletTrie(std::vector<cpp_int>::iterator row_begin, std::vector<cpp_int>::iterator row_end);

            template <class Container>
            WaveletTrie(Container&& rows);

            //destructor
            ~WaveletTrie() noexcept;

            cpp_int at(size_t i, size_t j = -1llu);

            size_t size() const;

            void insert(const WaveletTrie &wtr, size_t i = -1llu);
            void insert(WaveletTrie&& wtr, size_t i = -1llu);

            void print();

            template <typename T>
            void insert(const T &a, size_t i = -1llu);

            size_t serialize(std::ostream &out) const;
            size_t load(std::istream &in);

            bool operator==(const WaveletTrie &other) const;
            bool operator!=(const WaveletTrie &other) const;

        public:
            class Node;
        private:
            Node* root;
        //public:
        //    typedef std::pair<Node*, size_t> iterator;
    };

    class WaveletTrie::Node {
        friend class WaveletTrie;
        public:
            //empty constructor
            Node() { }

            Node(const size_t count);

            Node(const cpp_int &alpha, const size_t count);

            //destructor
            ~Node() noexcept;

            //Copy constructor
            Node(const Node &that);

            //move constructor
            Node(Node&& that) noexcept;

            Node& operator=(const Node &that);
            Node& operator=(Node&& that) noexcept;

            //Constructor from array
            //template <class Iterator>
            //Node(const Iterator &row_begin, const Iterator &row_end,
            //        const size_t &col = 0, Prefix prefix = Prefix());
            template <class Iterator>
            void fill_beta(const Iterator &row_begin, const Iterator &row_end,
                    const size_t &col, utils::ThreadPool &thread_queue, Prefix prefix = Prefix());

            size_t serialize(std::ostream &out) const;
            size_t load(std::istream &in);

            bool operator==(const Node &other) const;
            bool operator!=(const Node &other) const;

            size_t size() const { return beta_.size(); }

            void fill_left(bool rightside = false);

            bool check(bool ind);

            size_t rank1(const size_t i);

            size_t rank0(const size_t i);

            void print(std::ostream &out = std::cout) const;

        protected:
            alpha_t alpha_ = 1;
            beta_t beta_;
            rank1_t rank1_;
            Node *child_[2] = {NULL, NULL};
            size_t popcount = 0;
            bool support = false;

        private:
            static void merge_(Node *curnode, Node *othnode, size_t i, utils::ThreadPool &thread_queue);

            template <class IndexContainer>
            static size_t next_different_bit_(const IndexContainer &a, const IndexContainer &b,
                    size_t col = 0, size_t next_col = -1llu);

            static size_t next_different_bit_alpha(const Node &curnode, const Node &othnode);

            template <class Iterator>
            static Prefix longest_common_prefix(
                    const Iterator &row_begin, const Iterator &row_end, const size_t &col);

            //template <class Vector>
            //void set_beta_(const Vector &bv);

            template <class IndexContainer>
            void set_alpha_(const IndexContainer &indices, size_t col, size_t col_end);

            template <class IndexContainer>
            void set_alpha_(const IndexContainer &indices, size_t col);

            int move_label_down_(size_t length);

            static bool overlap_prefix_(Node &curnode, Node &othnode);

            static void merge_beta_(Node &curnode, const Node &othnode, size_t i = -1llu);

            void fill_ancestors(Node&& othnode, bool ind, const size_t i);

            bool is_leaf() const;

            template <class Container>
            static void push_child(
                    Container &nodes,
                    Node *curnode, Node *othnode,
                    bool ind, const size_t i
#ifndef NPRINT
                    , std::string path
#endif
            );

    };

};

#endif // __WAVELET_TRIE___
