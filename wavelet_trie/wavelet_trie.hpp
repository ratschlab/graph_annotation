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
#include "sdsl_utils.hpp"
#include "cpp_utils.hpp"


namespace annotate {

typedef cpp_int alpha_t;
typedef rrr_t beta_t;
//typedef bv_t beta_t;
typedef beta_t::rank_1_type rank1_t;

class Prefix {
  public:
    Prefix() { }
    Prefix(size_t col, bool allequal) : col(col), allequal(allequal) { }
    size_t col = -1llu;
    bool allequal = true;
};

class WaveletTrie {
  public:
    class Node;
    WaveletTrie();
    WaveletTrie(size_t p);

  private:
    WaveletTrie(const Node &node, size_t p = 1);
    WaveletTrie(Node *node, size_t p = 1);

  public:
    // copy constructor
    WaveletTrie(const WaveletTrie &other);

    // move constructor
    WaveletTrie(WaveletTrie&& other) noexcept;

    // assign
    WaveletTrie& operator=(const WaveletTrie& other);
    WaveletTrie& operator=(WaveletTrie&& other) noexcept;

    template <class Iterator>
    WaveletTrie(Iterator row_begin, Iterator row_end, size_t p = 1);
    //WaveletTrie(std::vector<cpp_int>::iterator row_begin, std::vector<cpp_int>::iterator row_end);

    template <class Container>
    WaveletTrie(Container&& rows, size_t p = 1);

    //destructor
    ~WaveletTrie() noexcept;

    cpp_int at(size_t i, size_t j = -1llu) const;

    void set_bit(size_t i, size_t j);

    void toggle_bit(size_t i, size_t j);

    void unset_bit(size_t i, size_t j);

    size_t size() const;

    void insert(const WaveletTrie &wtr, size_t i = -1llu);
    void insert(WaveletTrie&& wtr, size_t i = -1llu);

    void remove(size_t j);

    void print();

    template <typename T>
    void insert(const T &a, size_t i = -1llu);

    size_t serialize(std::ostream &out) const;
    size_t load(std::istream &in);

    bool operator==(const WaveletTrie &other) const;
    bool operator!=(const WaveletTrie &other) const;

    size_t get_p() const { return p_; }

    void set_p(size_t p) { p_ = p; }

  private:
    Node* root = NULL;
    size_t p_; // number of threads

    Node* traverse_down(Node *node, size_t &i, size_t &j);

}; // WaveletTrie

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

    void fill_left(size_t i, bool rightside);

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

    void fill_ancestors(Node &othnode, bool ind, const size_t i);

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

}; // WaveletTrie::Node

}; // annotate

#endif // __WAVELET_TRIE___
