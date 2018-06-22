#include "wavelet_trie.hpp"
#include <omp.h>
#include <thread>
#include <mutex>
#include <future>

namespace annotate {

std::mutex construct_mtx, merge_mtx;

WaveletTrie::WaveletTrie() : root(nullptr), p_(1) {}
WaveletTrie::WaveletTrie(size_t p) : root(nullptr), p_(p) {}
//WaveletTrie::WaveletTrie() : root(nullptr), thread_queue_(1) {}
//WaveletTrie::WaveletTrie(size_t p) : root(nullptr), thread_queue_(p) {}

//copy
WaveletTrie::WaveletTrie(const Node &node, size_t p)
    : root(new Node(const_cast<const Node&>(node))),
      //thread_queue_(p) {}
      p_(p) {}

WaveletTrie::WaveletTrie(Node *node, size_t p)
    : root(node),
      //thread_queue_(p) {}
      p_(p) {}

// copy constructor
WaveletTrie::WaveletTrie(const WaveletTrie &other)
    //: thread_queue_(other.thread_queue_.num_workers()) {
    : p_(other.p_) {
    if (other.root) {
        root = new Node(const_cast<const Node&>(*other.root));
    }
}

// move constructor
WaveletTrie::WaveletTrie(WaveletTrie&& other) noexcept
    : root(other.root),
      //thread_queue_(std::move(other.thread_queue_)) {
      p_(other.p_) {
    other.root = nullptr;
}

// copy assign
WaveletTrie& WaveletTrie::operator=(const WaveletTrie &other) {
    if (&other != this) {
        WaveletTrie tmp(const_cast<const WaveletTrie&>(other));
        std::swap(*this, tmp);
    }
    return *this;
}

// move assign
WaveletTrie& WaveletTrie::operator=(WaveletTrie&& other) noexcept {
    std::swap(root, other.root);
    std::swap(p_, other.p_);
    return *this;
}

template <class Iterator>
WaveletTrie::WaveletTrie(Iterator row_begin, Iterator row_end, size_t p)
    //: thread_queue_(p) {
    : p_(p) {
    if (std::distance(row_begin, row_end) > 0) {
        Prefix prefix = WaveletTrie::Node::longest_common_prefix(row_begin, row_end, 0);
        if (prefix.allequal) {
            root = new Node(std::distance(row_begin, row_end));
            root->set_alpha_(*row_begin, 0);
        } else {
            utils::ThreadPool thread_queue(p_);
            root = new Node();
            root->set_alpha_(*row_begin, 0, prefix.col);
            thread_queue.enqueue([=, &thread_queue]() {
                root->fill_beta(row_begin, row_end, 0, thread_queue, prefix);
            });
            thread_queue.join();
        }
    } else {
        root = NULL;
    }
#ifdef PRINT
    print();
    std::cout << "\n";
#endif
}


size_t WaveletTrie::serialize(std::ostream &out) const {
    if (root)
        return root->serialize(out);
    else
        return Node().serialize(out);
}

size_t WaveletTrie::load(std::istream &in) {
    if (root) {
        delete root;
    }
    root = new Node();
    size_t size = root->load(in);
    if (!root->beta_.size()) {
        delete root;
        root = NULL;
        return 0;
    }
    return size;
}

bool WaveletTrie::operator==(const WaveletTrie &other) const {
    if (size() != other.size())
        return false;
    if (static_cast<bool>(root) != static_cast<bool>(other.root))
        return false;
    if (root == nullptr)
        return true;
    return *root == *other.root;
}

bool WaveletTrie::operator!=(const WaveletTrie &other) const {
    return !(*this == other);
}

size_t WaveletTrie::Node::serialize(std::ostream &out) const {
    size_t written_bytes = ::annotate::serialize(out, alpha_)
                         + beta_.serialize(out);

    char val = (bool)child_[0] | ((uint8_t)((bool)child_[1]) << 1);
    out.write(&val, 1);
    written_bytes++;
    if (val == 1) {
        std::cerr << "ERROR: weird case\n";
        exit(1);
    }
    if (child_[0]) {
        assert(val & 1);
        written_bytes += child_[0]->serialize(out);
    }
    if (child_[1]) {
        assert(val & 2);
        written_bytes += child_[1]->serialize(out);
    }
    return written_bytes;
}

void WaveletTrie::Node::print(std::ostream &out) const {
    out << alpha_ << ":" << beta_ << ";" << is_leaf() << std::endl;
}

std::pair<size_t, size_t> WaveletTrie::stats() const {
    std::pair<size_t, size_t> num_uniq_set_bits = {0, 0};
    if (!root)
        return num_uniq_set_bits;
    std::stack<Node*> node_stack;
    node_stack.emplace(root);
    while (node_stack.size()) {
        Node *curnode = node_stack.top();
        node_stack.pop();
        num_uniq_set_bits.second += (popcount(curnode->alpha_) - 1) * curnode->size() + curnode->popcount;
        if (curnode->is_leaf()) {
            num_uniq_set_bits.first++;
            continue;
        }
        if (curnode->child_[0])
            node_stack.emplace(curnode->child_[0]);
        if (curnode->child_[1])
            node_stack.emplace(curnode->child_[1]);
    }
    return num_uniq_set_bits;
}

size_t WaveletTrie::Node::load(std::istream &in) {
    if (child_[0])
        delete child_[0];
    if (child_[1])
        delete child_[1];
    alpha_ = ::annotate::load(in);
    beta_.load(in);
    sdsl::util::init_support(rank1_, &beta_);
    support = true;
    popcount = beta_.size() ? rank1_(size()) : 0;
    char val;
    in.read(&val, 1);
    if (val >= '0') {
        val -= '0';
    }
    if (val == 1) {
        std::cerr << "ERROR: weird case\n";
        exit(1);
    }
    assert(!popcount == !val);
    if (val & 1) {
        //left child exists
        child_[0] = new Node();
        child_[0]->load(in);
    }
    if (val & 2) {
        //right child exists
        child_[1] = new Node();
        child_[1]->load(in);
    }
    return 0;
}

bool WaveletTrie::Node::operator==(const WaveletTrie::Node &other) const {
    if (alpha_ != other.alpha_) {
#ifdef PRINT
        print(); other.print();
#endif
        return false;
    }
    if (sdsl::util::to_string(beta_) != sdsl::util::to_string(other.beta_)) {
#ifdef PRINT
        print(); other.print();
#endif
        return false;
    }
    if (popcount != other.popcount) {
#ifdef PRINT
        print(); other.print();
#endif
        return false;
    }
    if ((bool)child_[0] != (bool)other.child_[0]) {
#ifdef PRINT
        std::cout << "left failed\n";
#endif
        return false;
    }
    if ((bool)child_[1] != (bool)other.child_[1]) {
#ifdef PRINT
        std::cout << "right failed\n";
#endif
        return false;
    }
    bool left_val = true;
    if (child_[0]) {
#ifdef PRINT
        std::cout << "left\n";
#endif
        left_val = *child_[0] == *other.child_[0];
    }
    if (!left_val)
        return false;
    if (child_[1]) {
#ifdef PRINT
        std::cout << "right\n";
#endif
        return *child_[1] == *other.child_[1];
    }
    return true;
}

bool WaveletTrie::Node::operator!=(const WaveletTrie::Node &other) const {
    return !(*this == other);
}

void WaveletTrie::print() const {
    if (!root) {
        std::cout << ":\t1:;1\n";
        return;
    }
    std::stack<std::pair<Node*,std::string>> nodestack;
    nodestack.emplace(root,"");
    while (nodestack.size()) {
        Node *curnode = nodestack.top().first;
        std::string parent = nodestack.top().second;
        nodestack.pop();
        std::cout << parent << ":\t" << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->is_leaf() << "\n";
        if (curnode->child_[0])
            nodestack.emplace(curnode->child_[0], parent + std::string("L"));
        if (curnode->child_[1])
            nodestack.emplace(curnode->child_[1], parent + std::string("R"));
    }
}


template<>
void WaveletTrie::Node::set_alpha_<cpp_int>(const cpp_int &alpha, pos_t col, pos_t col_end) {
    //alpha_ = (*row_begin >> col) % mask; //reference
    mpz_t& alph = alpha_.backend().data();
    mpz_tdiv_q_2exp(alph, alpha.backend().data(), col);
    clear_after(alph, col_end - col);
    bit_set(alph, col_end - col);
}

template <>
void WaveletTrie::Node::set_alpha_<std::set<pos_t>>(const std::set<pos_t> &indices, pos_t col, pos_t col_end) {
    alpha_ = 0;
    auto &alph = alpha_.backend().data();
    for (auto it = indices.lower_bound(col); it != indices.end(); ++it) {
        if (*it >= col_end)
            break;
        bit_set(alph, *it - col);
    }
    bit_set(alph, col_end - col);
}

template <class IndexContainer>
void WaveletTrie::Node::set_alpha_(const IndexContainer &indices, pos_t col, pos_t col_end) {
    alpha_ = 0;
    auto &alph = alpha_.backend().data();
    for (auto it = indices.begin(); it != indices.end(); ++it) {
        //if (*it >= col_end)
        //    break;
        if (*it >= col && *it < col_end)
            bit_set(alph, *it - col);
    }
    bit_set(alph, col_end - col);
}

template<>
void WaveletTrie::Node::set_alpha_<cpp_int>(const cpp_int &alpha, pos_t col) {
    mpz_t& alph = alpha_.backend().data();
    mpz_tdiv_q_2exp(alph, alpha.backend().data(), col);
    if (is_nonzero(alpha_)) {
        bit_set(alpha_, msb(alpha_) + 1);
    } else {
        bit_set(alpha_, 0);
    }
}

template <>
void WaveletTrie::Node::set_alpha_<std::set<pos_t>>(const std::set<pos_t> &indices, pos_t col) {
    alpha_ = 0;
    auto &alph = alpha_.backend().data();
    for (auto it = indices.lower_bound(col); it != indices.end(); ++it) {
        bit_set(alph, *it - col);
    }
    if (is_nonzero(alpha_)) {
        bit_set(alpha_, msb(alpha_) + 1);
    } else {
        bit_set(alpha_, 0);
    }
}

template <class IndexContainer>
void WaveletTrie::Node::set_alpha_(const IndexContainer &indices, pos_t col) {
    alpha_ = 0;
    auto &alph = alpha_.backend().data();
    for (auto it = indices.begin(); it != indices.end(); ++it) {
        if (*it >= col)
            bit_set(alph, *it - col);
    }
    if (is_nonzero(alpha_)) {
        bit_set(alpha_, msb(alpha_) + 1);
    } else {
        bit_set(alpha_, 0);
    }
}

template <class Iterator>
void WaveletTrie::Node::fill_beta(const Iterator &row_begin, const Iterator &row_end,
        const pos_t &col, utils::ThreadPool &thread_queue, Prefix prefix) {
    //TODO col already used by Prefix?
    std::ignore = col;
    if (std::distance(row_begin, row_end)) {
        assert(prefix.col != static_cast<pos_t>(-1));
        assert(!prefix.allequal);
        pos_t col_end = prefix.col;

        //set beta and compute common prefices
        bv_t beta;
        beta.resize(std::distance(row_begin, row_end));
        Prefix prefices[2];
        Iterator split = row_begin;
        std::vector<typename std::iterator_traits<Iterator>::value_type> right_children;
        sdsl::util::set_to_value(beta, 0);
        auto begin = std::make_move_iterator(row_begin);
        auto end = std::make_move_iterator(row_end);
        for (auto it = begin; it != end; ++it) {
            if (bit_test(*it, col_end)) {
                beta[std::distance(begin, it)] = 1;
                right_children.emplace_back(*it);
                prefices[1].col = next_different_bit_(
                        right_children.front(), right_children.back(),
                        col_end + 1, prefices[1].col);
            } else {
                if (std::distance(row_begin, split) != std::distance(begin, it)) //prevent setting equality on same object
                    *split = *it;
                prefices[0].col = next_different_bit_(
                        *row_begin, *split,
                        col_end + 1, prefices[0].col);
                split++;
            }
        }
        popcount = right_children.size();
        beta_ = beta_t(beta);
        support = false;
        //assert(popcount == rank1(size()));

        //distribute to left and right children
        assert(split != row_begin && split != row_end);
        std::move(right_children.begin(), right_children.end(), split);
        right_children.clear();

        //TODO: copied code here
        //handle trivial cases first
        if (prefices[0].col == static_cast<pos_t>(-1)) {
            child_[0] = new Node(std::distance(row_begin, split));
            child_[0]->set_alpha_(*row_begin, col_end + 1);
            assert(child_[0]->size() == rank0(beta_.size()));
        }

        if (prefices[1].col == static_cast<pos_t>(-1)) {
            child_[1] = new Node(std::distance(split, row_end));
            child_[1]->set_alpha_(*split, col_end + 1);
            assert(child_[1]->size() == rank1(beta_.size()));
        }

        //then recursive calls
        if (prefices[0].col != static_cast<pos_t>(-1)) {
            prefices[0].allequal = false;
            child_[0] = new Node();
            child_[0]->set_alpha_(*row_begin, col_end + 1, prefices[0].col);
            std::lock_guard<std::mutex> lock(construct_mtx);
            thread_queue.enqueue([=, &thread_queue]() {
                child_[0]->fill_beta(row_begin, split, col_end + 1, thread_queue, prefices[0]);
            });
        }

        if (prefices[1].col != static_cast<pos_t>(-1)) {
            prefices[1].allequal = false;
            child_[1] = new Node();
            child_[1]->set_alpha_(*split, col_end + 1, prefices[1].col);
            std::lock_guard<std::mutex> lock(construct_mtx);
            thread_queue.enqueue([=, &thread_queue]() {
                child_[1]->fill_beta(split, row_end, col_end + 1, thread_queue, prefices[1]);
            });
        }
    }
}

WaveletTrie::Node::Node(const size_t count)
  : beta_(beta_t(bv_t(count))), popcount(0), support(false) {}

WaveletTrie::Node::Node(const cpp_int &alpha, const size_t count)
  : Node(count) {
    alpha_ = alpha;
}

WaveletTrie::Node& WaveletTrie::Node::operator=(const Node &that) {
    if (&that != this) {
        Node tmp(that);
        std::swap(*this, tmp);
    }
    return *this;
}

WaveletTrie::Node& WaveletTrie::Node::operator=(Node&& that) noexcept {
    std::swap(alpha_, that.alpha_);
    std::swap(beta_, that.beta_);
    std::swap(rank1_, that.rank1_);
    std::swap(popcount, that.popcount);
    std::swap(support, that.support);
    std::swap(child_[0], that.child_[0]);
    std::swap(child_[1], that.child_[1]);
    that.child_[0] = nullptr;
    that.child_[1] = nullptr;
    return *this;
}

WaveletTrie::Node::Node(const Node &that)
    : alpha_(that.alpha_), beta_(that.beta_),
      rank1_(that.rank1_),
      popcount(that.popcount),
      support(that.support) {
    if (that.child_[0]) {
        child_[0] = new Node(const_cast<const Node&>(*that.child_[0]));
    }
    if (that.child_[1]) {
        child_[1] = new Node(const_cast<const Node&>(*that.child_[1]));
    }
}

WaveletTrie::Node::Node(Node&& that) noexcept
    : alpha_(std::move(that.alpha_)), beta_(std::move(that.beta_)),
      rank1_(std::move(that.rank1_)),
      popcount(that.popcount),
      support(that.support) {
    child_[0] = that.child_[0];
    child_[1] = that.child_[1];
    that.child_[0] = nullptr;
    that.child_[1] = nullptr;
}

template <class Container>
WaveletTrie::WaveletTrie(Container&& rows, size_t p) : WaveletTrie::WaveletTrie(rows.begin(), rows.end(), p) {}

WaveletTrie::~WaveletTrie() noexcept {
    delete root;
}

//TODO: use callbacks to avoid duplicated code
WaveletTrie::Node* WaveletTrie::traverse_down(Node *node, size_t &i, pos_t &j) {
    size_t curmsb;
    while (!node->is_leaf() && (curmsb = msb(node->alpha_) + 1) < j) {
        j -= curmsb;
        if (node->beta_[i]) {
            assert(node->child_[1]);
            i = node->rank1(i);
            node = node->child_[1];
        } else {
            assert(node->child_[0]);
            i = node->rank0(i);
            node = node->child_[0];
        }
    }
    return node;
}

void WaveletTrie::toggle_bit(size_t i, pos_t j) {
    assert(i < size());
    WaveletTrie wtr_int(traverse_down(root, i, j));
    auto cur_annot = wtr_int.at(i);
    if (bit_test(cur_annot, j)) {
        bit_unset(cur_annot, j);
    } else {
        bit_set(cur_annot, j);
    }
    wtr_int.insert(WaveletTrie(std::vector<cpp_int>{cur_annot}), i);
    wtr_int.remove(i + 1);
    wtr_int.root = NULL;
}

cpp_int WaveletTrie::at(size_t i, pos_t j) const {
    assert(i < size());
    Node *node = root;
    size_t length = 0;
    cpp_int annot;
    if (!node)
        return annot;
    while (!node->is_leaf() && length < j) {
        annot |= node->alpha_ << length;
        length += msb(node->alpha_) + 1;
        if (node->beta_[i]) {
            assert(node->child_[1]);
            i = node->rank1(i);
            node = node->child_[1];
        } else {
            bit_unset(annot, length - 1);
            assert(node->child_[0]);
            i = node->rank0(i);
            node = node->child_[0];
        }
    }
    annot |= node->alpha_ << length;
    bit_unset(annot, msb(annot));
    return annot;
}

void WaveletTrie::set_bit(size_t i, pos_t j) {
    assert(i < size());
    WaveletTrie wtr_int(traverse_down(root, i, j), p_);
    auto cur_annot = wtr_int.at(i);
    if (bit_test(cur_annot, j)) {
        wtr_int.root = NULL;
        return;
    }
    bit_set(cur_annot, j);

    WaveletTrie wtr_int_ins(std::vector<cpp_int>{cur_annot});
    wtr_int.insert(std::move(wtr_int_ins), i);
    wtr_int.remove(i + 1);
    wtr_int.root = NULL;
}

template <class Container>
void WaveletTrie::set_bits(Container &is, pos_t j) {
    if (is.empty())
        return;
    assert(*is.rbegin() < size());

    std::vector<cpp_int> edges_old;
    std::vector<pos_t> from, to;
    from.reserve(is.size());
    to.reserve(is.size());
    size_t oldsize = size();
    for (auto &i : is) {
        from.push_back(size() + edges_old.size());
        to.push_back(i);
        edges_old.push_back(at(i));
        bit_set(edges_old.back(), j);
    }
    insert(WaveletTrie(edges_old));
    swap(from, to, oldsize);
    //std::cout << "Swap\n";
    std::vector<pos_t> indices(edges_old.size());
    std::iota(indices.begin(), indices.end(), oldsize);
    remove(indices);
    assert(size() == oldsize);
}
template void WaveletTrie::set_bits(std::vector<pos_t>&, pos_t);
template void WaveletTrie::set_bits(std::set<pos_t>&, pos_t);

void WaveletTrie::unset_bit(size_t i, pos_t j) {
    assert(i < size());
    WaveletTrie wtr_int(traverse_down(root, i, j));
    auto cur_annot = wtr_int.at(i);
    if (!bit_test(cur_annot, j)) {
        wtr_int.root = NULL;
        return;
    }
    bit_unset(cur_annot, j);

    WaveletTrie wtr_int_ins(std::vector<cpp_int>{cur_annot});
    wtr_int.insert(std::move(wtr_int_ins), i);
    wtr_int.remove(i + 1);
    wtr_int.root = NULL;
}

size_t WaveletTrie::size() const {
    return root ? root->size() : 0;
}

bool WaveletTrie::Node::is_leaf() const {
    return (popcount == 0);
}

bool WaveletTrie::Node::check(bool ind) {
    size_t rank = ind ? popcount : size() - popcount;
    //assert(rank1(size()) == popcount);
    //assert(popcount != size());
    if (!popcount) {
        assert(!child_[1]);
        assert(!child_[0]);
    }
    if (child_[ind]) {
        if (is_leaf())
            return false;
        assert(child_[ind]->size() == rank);
        if (child_[ind]->size() != rank)
            return false;
        return child_[ind]->check(0) && child_[ind]->check(1);
    } else {
        assert(is_leaf() || (child_[!ind] && rank == 0));
        if (!is_leaf() && (!child_[!ind] || rank > 0))
            return false;
    }
    return true;
}

void WaveletTrie::Node::fill_left(size_t i, bool rightside) {
    size_t lrank = size() - popcount;
    assert(lrank == rank0(size()));
    Node *jnode = this;
    if (lrank) {
        while (jnode->child_[0]) {
            //assert(lrank);
            Node *lchild = jnode->child_[0];
            assert(lrank >= lchild->size());
            if (lrank == lchild->size())
                return;
            lchild->move_label_down_(lsb(lchild->alpha_));
            assert(lchild->alpha_ == 1 || ((lchild->alpha_ | 1) == lchild->alpha_ + 1));
            if (rightside) {
#ifdef PRINT
                std::cout << "bar\t" << i << "\n";
#endif
                assert(i <= lchild->size());
#ifdef PRINT
                std::cout << lchild->beta_ << "\n";
#endif
                lchild->beta_ = beta_t(insert_zeros(lchild->beta_,
                                                lrank - lchild->beta_.size(),
                                                i));
                lchild->support = false;
                i = lchild->rank0(i);
#ifdef PRINT
                std::cout << lchild->beta_ << "\n";
#endif
            } else {
#ifdef PRINT
                std::cout << "foo\t" << i << "\n";
#endif
                assert(lrank >= i + lchild->size());
                auto bv = insert_zeros(lchild->beta_,
                        i,
                        0);
                lchild->beta_ = beta_t(insert_zeros(bv,
                            lrank - i - lchild->size(),
                            bv.size()));
                lchild->support = false;
            }
            lrank -= lchild->popcount;
            jnode = jnode->child_[0];
        }
        assert(lrank);
        assert(jnode->popcount || (!jnode->popcount && jnode->size() == lrank));
        if (jnode->popcount) {
            assert(lrank + jnode->popcount == jnode->size());
            jnode->child_[0] = new Node(lrank);
        }
    }
}

void WaveletTrie::Node::fill_ancestors(Node &othnode, bool ind, const size_t i) {
    if (child_[ind]) {
        if (!ind
            && othnode.is_leaf()
            && !othnode.child_[0]
            && !othnode.child_[1]
            && (popcount - othnode.popcount)) {
            assert(popcount);
            assert(othnode.popcount == 0);
#ifdef PRINT
            std::cout << (ind ? "R" : "L") << i << "\n";
#endif
            fill_left(i, true);
        }
    } else {
        assert(child_[!ind] || (popcount == othnode.popcount));
        if (othnode.child_[ind]) {
            std::swap(child_[ind], othnode.child_[ind]);
            if (!ind
                && (popcount == othnode.popcount)
                && !othnode.is_leaf()) {
#ifdef PRINT
                std::cout << (ind ? "R" : "L") << "s" << i << "\n";
#endif
                fill_left(i, false);
            }
            assert(popcount);
            assert(popcount == rank1(size()));
        } else if (!ind && popcount && size() > popcount) {
            assert(size() - popcount == rank0(size()));
            child_[ind] = new Node(size() - popcount);
        }
    }
}

void WaveletTrie::Node::merge_(Node *curnode, Node *othnode, size_t i, utils::ThreadPool &thread_queue) {
    assert(curnode);
    assert(othnode);
    assert(curnode->size());
    assert(othnode->size());
    assert(i <= curnode->size());

    while (curnode && othnode) {
        assert(curnode->size());
        assert(i <= curnode->size());
        assert(curnode->check(0));
        assert(curnode->check(1));
#ifdef PRINT
        std::cout << "" << "\t" << i << "\t"
                  << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->is_leaf() << "\t"
                  << othnode->alpha_ << ":" << othnode->beta_ << ";" << othnode->is_leaf() << "\t->\t";
#endif
        Node::overlap_prefix_(*curnode, *othnode);

        //update insertion point and merge betas
        assert(othnode);
        Node::merge_beta_(*curnode, *othnode, i);
#ifdef PRINT
        std::cout << curnode->alpha_ << ":" << curnode->beta_ << ";" << curnode->is_leaf() << "\t"
                  << othnode->alpha_ << ":" << othnode->beta_ << ";" << othnode->is_leaf() << "\n";
#endif
        bool left  = curnode->child_[0] && othnode->child_[0];
        bool right = curnode->child_[1] && othnode->child_[1];
#ifdef PRINT
        std::cout << "P" << i << "\n";
#endif
        size_t il = curnode->rank0(i);
        size_t ir = curnode->rank1(i);
        curnode->fill_ancestors(*othnode, 0, il);
        curnode->fill_ancestors(*othnode, 1, ir);
        if ((bool)curnode->child_[0] != (bool)curnode->child_[1]) {
            std::cerr << "ERROR: merging imbalanced error" << std::endl;
            exit(1);
        }

        if (left && right) {
            //send smaller problem to own thread
            if (std::min(curnode->child_[0]->size(),
                         othnode->child_[0]->size())
              > std::min(curnode->child_[1]->size(),
                         othnode->child_[1]->size())) {
                std::lock_guard<std::mutex> lock(merge_mtx);
                thread_queue.enqueue([=, &thread_queue]() {
                    merge_(curnode->child_[1], othnode->child_[1], ir, thread_queue);
                });
                curnode = curnode->child_[0];
                othnode = othnode->child_[0];
                i = il;
            } else {
                std::lock_guard<std::mutex> lock(merge_mtx);
                thread_queue.enqueue([=, &thread_queue]() {
                    merge_(curnode->child_[0], othnode->child_[0], il, thread_queue);
                });
                curnode = curnode->child_[1];
                othnode = othnode->child_[1];
                i = ir;
            }
        } else if (left) {
            curnode = curnode->child_[0];
            othnode = othnode->child_[0];
            i = il;
        } else if (right) {
            curnode = curnode->child_[1];
            othnode = othnode->child_[1];
            i = ir;
        } else {
            curnode = NULL;
            othnode = NULL;
        }
    }
}

void WaveletTrie::insert(const WaveletTrie &wtr, size_t i) {
    if (!wtr.root) {
        return;
    }
    if (!root) {
        root = new Node(const_cast<const Node&>(*wtr.root));
        return;
    }
    if (i == static_cast<size_t>(-1)) {
        i = size();
    }
    WaveletTrie tmp(wtr);
    utils::ThreadPool thread_queue(p_);
    thread_queue.enqueue([=, &thread_queue, &tmp]() {
        Node::merge_(root, tmp.root, i, thread_queue);
    });
    thread_queue.join();
}

void WaveletTrie::insert(WaveletTrie&& wtr, size_t i) {
    if (!wtr.root) {
        return;
    }
    if (!root) {
        std::swap(root, wtr.root);
        return;
    }
    if (i == -1llu) {
        i = size();
    }
    utils::ThreadPool thread_queue(p_);
    thread_queue.enqueue([=, &thread_queue, &wtr]() {
        Node::merge_(root, wtr.root, i, thread_queue);
    });
    thread_queue.join();
}

template <typename T>
void WaveletTrie::insert(const T &a, size_t i) {
    Node *next = new Node(&a, &a + 1, 0);
    if (!i) {
        root = next;
    } else {
        WaveletTrie wtr;
        wtr.root = next;
        insert(wtr, i);
    }
}

void WaveletTrie::swap(const std::vector<pos_t> &from,
                       const std::vector<pos_t> &to,
                       pos_t remove_after) {
    if (from.size() != to.size()) {
        std::cerr << "Invalid map\n";
        exit(1);
    }
    if (!from.size())
        return;
    std::ignore = remove_after;

    //std::cout << "Base\n";
    //print();
    //std::cout << "\n";
    typedef std::pair<pos_t, pos_t> SwapState;
    struct NodeState {
        Node *node;
        std::vector<SwapState> swap_states;
        std::vector<SwapState> move_states;
        size_t remove_after;
        NodeState(Node *node) : node(node) {}
        NodeState(Node *node,
                  const std::vector<SwapState> &swap_states,
                  const std::vector<SwapState> &move_states = {})
            : node(node), swap_states(swap_states), move_states(move_states) {}
    };

    utils::ThreadPool thread_queue(p_);
    std::stack<NodeState> node_stack;
    std::vector<SwapState> init_states;
    init_states.reserve(from.size());
    for (size_t i = 0; i < from.size(); ++i) {
        init_states.emplace_back(from[i], to[i]);
    }
    node_stack.emplace(root, init_states);
    //node_stack.top().remove_after = remove_after;

    while (node_stack.size()) {
        auto curstate = node_stack.top();
        //auto curstate = node_stack.top().get();
        node_stack.pop();
        Node *node = curstate.node;
        assert(node);
        if (node->is_leaf() || (curstate.swap_states.empty() && curstate.move_states.empty()))
            continue;

        NodeState next_states[2] = {NodeState{node->child_[0]}, NodeState{node->child_[1]}};
        //next_states[1].remove_after = node->rank1(curstate.remove_after);
        //next_states[0].remove_after = curstate.remove_after - next_states[1].remove_after;
        for (auto &state : curstate.swap_states) {
            assert(state.first < node->size());
            assert(state.second < node->size());
            bool from_child = node->beta_[state.first];
            bool to_child = node->beta_[state.second];
            if (from_child == to_child) {
                //stay in swap state
                if (!node->child_[from_child]->is_leaf()) {
                    next_states[from_child].swap_states.emplace_back(
                            from_child ? node->rank1(state.first)  : node->rank0(state.first),
                            from_child ? node->rank1(state.second) : node->rank0(state.second));
                }
            } else {
                //transfer to move state
                if (!node->child_[from_child]->is_leaf()) {
                   next_states[from_child].move_states.emplace_back(
                           from_child ? node->rank1(state.first)  : node->rank0(state.first),
                           from_child ? node->rank1(state.second) : node->rank0(state.second));
                }
                if (!node->child_[to_child]->is_leaf()) {
                   next_states[to_child].move_states.emplace_back(
                           to_child ? node->rank1(state.second): node->rank0(state.second),
                           to_child ? node->rank1(state.first) : node->rank0(state.first));
                }
            }
        }
        for (auto &state : curstate.move_states) {
            bool from_child = node->beta_[state.first];
            if (!node->child_[from_child]->is_leaf()) {
                next_states[from_child].move_states.emplace_back(
                        from_child ? node->rank1(state.first)  : node->rank0(state.first),
                        from_child ? node->rank1(state.second) : node->rank0(state.second));
            }
        }
        node_stack.emplace(std::move(next_states[0]));
        node_stack.emplace(std::move(next_states[1]));

        thread_queue.enqueue([=]() {
            auto bv = swap_bits(node->beta_, curstate.swap_states);
            node->beta_ = beta_t(move_bits(bv,
                                           curstate.move_states,
                                           curstate.remove_after));
            node->support = false;
        });
    }
    thread_queue.join();
}

void WaveletTrie::remove(pos_t j) {
    remove(std::vector<pos_t>{j});
}

void WaveletTrie::remove(const std::vector<pos_t> &js) {
    if (!js.size())
        return;
    assert(js.back() < size());
#ifndef NDEBUG
    size_t expsize = size() - js.size();
#endif
    //if (size() == 1) {
    if (size() <= js.size()) {
        delete root;
        root = NULL;
        assert(size() == 0);
        return;
    }

    struct NodeState {
        Node *node;
        std::vector<pos_t> js;
        NodeState(Node *node, const std::vector<pos_t> &js)
            : node(node), js(std::move(js)) {}
    };

    std::stack<NodeState> node_stack;
    node_stack.emplace(root, std::move(js));

    //Node *node = root;
    //while (!node->is_leaf()) {
    while (node_stack.size()) {
        auto curstate = node_stack.top();
        node_stack.pop();
        Node *node = curstate.node;
        assert(node);
        //while (true) {
        if (node->is_leaf()) {
            assert(!node->popcount);
            node->beta_ = beta_t(bv_t(node->size() - curstate.js.size()));
            //curstate.node->beta_ = beta_t(bv_t(curstate.node->size() - std::distance(curstate.begin, curstate.end)));
            node->support = false;
            continue;
            //break;
        }
        //bool curbit = node->beta_[j];
        //size_t next_j;
        size_t curpopcount = node->popcount;
        std::vector<pos_t> next_js[2];
        //std::vector<size_t> js_temp(curstate.begin, curstate.end);
        //std::vector<size_t> next_js0;
        for (size_t i = 0; i < curstate.js.size(); ++i) {
            if (node->beta_[curstate.js[i]]) {
            //if (node->beta_[*it]) {
                next_js[1].emplace_back(node->rank1(curstate.js[i]));
                curpopcount--;
            } else {
                next_js[0].emplace_back(node->rank0(curstate.js[i]));
            }
        }
        //std::cout << "\n";
        /*
        if (curbit) {
            next_j = node->rank1(j);
            curpopcount--;
        } else {
            next_j = node->rank0(j);
        }
        */
        //if ((curbit && !curpopcount)
        //        || (!curbit && curpopcount == node->beta_.size() - 1)) {
        //curpopcount += added - subtracted;
        //fold one side in
        if (!curpopcount
                || curpopcount == node->beta_.size() - curstate.js.size()) {
            //merge nodes
            bool child = static_cast<bool>(curpopcount);
            delete node->child_[!child];
            node->child_[!child] = NULL;
            size_t curmsb = msb(node->alpha_);
            if (!child) {
                bit_unset(node->alpha_, curmsb);
            }
            Node *othnode = node->child_[child];
            assert(othnode);
            assert(othnode->size() >= next_js[child].size());
            node->child_[child] = NULL;
            node->alpha_ |= othnode->alpha_ << (curmsb + 1);
            assert(node->alpha_ != 0);
            node->beta_ = std::move(othnode->beta_);
            node->popcount = othnode->popcount;
            node->rank1_.set_vector(&node->beta_);
            std::swap(node->child_[0], othnode->child_[0]);
            std::swap(node->child_[1], othnode->child_[1]);
            delete othnode;
            if (node->is_leaf()) {
                curmsb = msb(node->alpha_);
                bit_unset(node->alpha_, curmsb);
                if (node->alpha_ != 0)
                    bit_set(node->alpha_, msb(node->alpha_) + 1);
                else
                    node->alpha_ = 1;
            }
            if (next_js[child].size()) {
                node_stack.emplace(node, std::move(next_js[child]));
                //curstate.js = std::move(next_js[child]);
            }
            //return;
            continue;
        } else {
            /*
            node->beta_ = beta_t(remove_range(node->beta_, j, j + 1));
            node->support = false;
            if (curbit)
                node->popcount--;
            */
            node->beta_ = beta_t(remove_bits(node->beta_, curstate.js));
            node->support = false;
            node->popcount = curpopcount;
        }
        //j = next_j;
        //node = curbit ? node->child_[1] : node->child_[0];
        //assert(node);
        /*
        if (next_js[1].size() && next_js[0].size()) {
            node_stack.emplace(node->child_[1], std::move(next_js[1]), std::move(next_swap[1]));
            node_stack.emplace(node->child_[0], std::move(next_js[0]), std::move(next_swap[0]));
            break;
        } else if (!next_js[1].size() && !next_js[0].size()) {
            break;
        } else {
            if (next_js[1].size()) {
                node = node->child_[1];
                curstate.js = std::move(next_js[1]);
                curstate.swap = std::move(next_swap[1]);
            }
            if (next_js[0].size()) {
                node = node->child_[0];
                curstate.js = std::move(next_js[0]);
                curstate.swap = std::move(next_swap[0]);
            }
        }
        */
        if (node->child_[1]->is_leaf()) {
            Node *child = node->child_[1];
            child->beta_ = beta_t(bv_t(curpopcount));
            child->support = false;
        } else if (next_js[1].size()) {
            node_stack.emplace(node->child_[1], std::move(next_js[1]));
        }
        if (node->child_[0]->is_leaf()) {
            Node *child = node->child_[0];
            child->beta_ = beta_t(bv_t(node->size() - curpopcount));
            child->support = false;
        } else if (next_js[0].size()) {
            node_stack.emplace(node->child_[0], std::move(next_js[0]));
        }
        //}
    }
    /*
    assert(!node->popcount);
    node->beta_ = beta_t(bv_t(node->size() - 1));
    node->support = false;
    */
    assert(expsize  == size());
}

WaveletTrie::Node::~Node() noexcept {
    if (child_[0]) {
        delete child_[0];
    }
    if (child_[1]) {
        delete child_[1];
    }
}

//return true if a change happened
bool WaveletTrie::Node::overlap_prefix_(Node &curnode, Node &othnode) {
    assert(curnode.alpha_ && othnode.alpha_);

    if (curnode.alpha_ == othnode.alpha_) {
        return false;
    }
    pos_t common_pref = next_different_bit_alpha(curnode, othnode);
#ifdef PRINT
    std::cout << common_pref << "\t";
#endif
    int cur = curnode.move_label_down_(common_pref);
#ifdef PRINT
    std::cout << curnode.alpha_ << ":" << curnode.beta_ << ";" << curnode.is_leaf();
    if (cur > -1) {
        std::cout << "," << curnode.child_[cur]->alpha_ << ":" << curnode.child_[cur]->beta_ << ";" << curnode.child_[cur]->is_leaf();
    }
    std::cout << "\t";
#endif
    int oth = othnode.move_label_down_(common_pref);
#ifdef PRINT
    std::cout << othnode.alpha_ << ":" << othnode.beta_ << ";" << othnode.is_leaf();
    if (oth > -1) {
        std::cout << "," << othnode.child_[oth]->alpha_ << ":" << othnode.child_[oth]->beta_ << ";" << othnode.child_[oth]->is_leaf();
    }
    std::cout << "\t->\t";
#endif
    assert(curnode.alpha_ == othnode.alpha_);
    if (cur > -1 && oth > -1 && curnode.child_[0] && othnode.child_[0]) {
        std::cerr << "ERROR: extra zero bit" << std::endl;
        exit(1);
    }
    if (cur > -1 && oth > -1 && curnode.child_[1] && othnode.child_[1]) {
        std::cerr << "ERROR: extra one bit" << std::endl;
        exit(1);
    }
    return true;
}

//find first different bit between two cpp_ints
template <class IndexContainer>
pos_t WaveletTrie::Node::next_different_bit_(
        const IndexContainer &a, const IndexContainer &b,
        pos_t col, pos_t next_col) {
    if (col == next_col)
        return next_col;
    pos_t ranges[3] = {next_bit(a, col), next_bit(b, col), next_col};
    while (ranges[0] == ranges[1] && ranges[0] < ranges[2] && ranges[1] < ranges[2]) {
        ranges[0] = next_bit(a, ranges[0] + 1);
        ranges[1] = next_bit(b, ranges[1] + 1);
    }
    pos_t next_col2 = std::min(ranges[0], ranges[1]);
    if (ranges[0] == ranges[1] || next_col2 >= ranges[2])
        next_col2 = ranges[2];
    return next_col2;
}

template<>
pos_t WaveletTrie::Node::next_different_bit_<cpp_int>(
        const cpp_int &a, const cpp_int &b,
        pos_t col, pos_t next_col) {
    if (col == next_col)
        return next_col;

    //try to move forward one limb at a time
    //get amount to shift based on limb size
    size_t shift = __builtin_clzll(mp_bits_per_limb) ^ 63;
    pos_t i = col >> shift;
    pos_t end = next_col >> shift;
    if (end > i) {
        const mpz_t &a_m = a.backend().data();
        const mpz_t &b_m = b.backend().data();
        auto a_l = mpz_limbs_read(a_m) + i;
        auto b_l = mpz_limbs_read(b_m) + i;
        auto a_e = mpz_limbs_read(a_m) + mpz_size(a_m);
        auto b_e = mpz_limbs_read(b_m) + mpz_size(b_m);
        if (a_l < a_e && b_l < b_e) {
            if (*a_l == *b_l) {
                col = (i + 1) << shift;
                ++a_l;
                ++b_l;
            } else if ((i << shift) < col) {
                //if col is in the middle of a limb, mask out earlier bits
                size_t mask = ~((1llu << (col % (1llu << shift))) - 1);
                if ((*a_l & mask) == (*b_l & mask)) {
                    ++a_l;
                    ++b_l;
                    col = (i + 1) << shift;
                }
            }
            while (a_l != a_e && b_l != b_e && *a_l == *b_l && i < end) {
                col += 1llu << shift;
                ++a_l;
                ++b_l;
                ++i;
            }
        }
    }

    //when at first different limb, use mpz API
    pos_t ranges[3] = {next_bit(a, col), next_bit(b, col), next_col};
    while (ranges[0] == ranges[1] && ranges[0] < ranges[2] && ranges[1] < ranges[2]) {
        ranges[0] = next_bit(a, ranges[0] + 1);
        ranges[1] = next_bit(b, ranges[1] + 1);
    }
    pos_t next_col2 = std::min(ranges[0], ranges[1]);
    if (ranges[0] == ranges[1] || next_col2 >= ranges[2])
        next_col2 = ranges[2];
    return next_col2;
}

pos_t WaveletTrie::Node::next_different_bit_alpha(const Node &curnode, const Node &othnode) {
    assert(curnode.alpha_ != 0);
    assert(othnode.alpha_ != 0);

    //TODO: replace this with a single pass
    auto cur = curnode.alpha_;
    auto oth = othnode.alpha_;
    pos_t curmsb = msb(cur);
    pos_t othmsb = msb(oth);
    bit_unset(cur, curmsb);
    bit_unset(oth, othmsb);
    pos_t next_set_bit = next_different_bit_(cur, oth);
    if (next_set_bit == static_cast<pos_t>(-1)) {
        if (curnode.is_leaf() == othnode.is_leaf()) {
            return std::min(curmsb, othmsb);
        }
        return std::max(curmsb, othmsb);
    } else {
        if (next_set_bit > curmsb && !curnode.is_leaf()) {
            next_set_bit = curmsb;
        }
        if (next_set_bit > othmsb && !othnode.is_leaf()) {
            next_set_bit = othmsb;
        }
    }
    return next_set_bit;
}

template <class Iterator>
Prefix WaveletTrie::Node::longest_common_prefix(const Iterator &row_begin, const Iterator &row_end, const pos_t &col) {
    Prefix prefix;
    //empty prefix
    if (row_begin >= row_end) {
        prefix.col = col;
        prefix.allequal = true;
        return prefix;
    }
    //prefix.col = static_cast<pos_t>(-1);
    for (auto it = row_begin + 1; it != row_end; ++it) {
        prefix.col = next_different_bit_(*row_begin, *it, col, prefix.col);
        if (prefix.col == col)
            break;
    }
    if (prefix.col == static_cast<pos_t>(-1)) {
        //all zeros or all equal
        if (is_nonzero(*row_begin)) {
            prefix.col = std::max(msb(*row_begin), col);
        } else {
            prefix.col = col;
        }
        return prefix;
    }
    prefix.allequal = false;
    return prefix;
}

int WaveletTrie::Node::move_label_down_(size_t length) {
    size_t len = msb(alpha_);
    if (length > len) {
        bit_unset(alpha_, len);
        bit_set(alpha_, length);
    } else if (length < len) {
        //split alpha and compute new beta
        //TODO: clean up
        Node *child = new Node();
        mpz_t& child_alpha = child->alpha_.backend().data();
        mpz_tdiv_q_2exp(child_alpha, alpha_.backend().data(), length + 1);
        //TODO: replace with operator=
        //child->set_beta_(beta_);
        child->beta_ = beta_;
        child->support = false;
        child->child_[0] = child_[0];
        child->child_[1] = child_[1];
        child->popcount = popcount;
        bool beta_bit = bit_test(alpha_, length);
        //set_beta_(bv_t(size(), beta_bit));
        beta_ = beta_t(bv_t(size(), beta_bit));
        support = false;
        //only want length bits left
        clear_after(alpha_, length);
        bit_set(alpha_, length);
        if (beta_bit) {
            child_[0] = NULL;
            child_[1] = child;
            popcount = size();
            assert(beta_.size() == child_[1]->size());
        } else {
            child_[0] = child;
            child_[1] = NULL;
            popcount = 0;
            assert(beta_.size() == child_[0]->size());
        }
        assert(popcount == rank1(size()));
        assert(msb(alpha_) == length);
        assert((child_[0] == NULL) ^ (child_[1] == NULL));
    }
    if (length >= len)
        return -1;
    if (child_[0])
        return 0;
    return 1;
}

void WaveletTrie::Node::merge_beta_(Node &curnode, const Node &othnode, size_t i) {
    assert(curnode.alpha_ == othnode.alpha_);
    assert(curnode.popcount == curnode.rank1(curnode.size()));
    if (!othnode.size()) {
        return;
    }
    // TODO: this test doesn't compile since othnode is now const
    //assert(othnode.popcount == othnode.rank1(othnode.size()));
    if (i == static_cast<size_t>(-1)) {
        i = curnode.beta_.size();
    }
    assert(i <= curnode.beta_.size());
    curnode.popcount += othnode.popcount;
    //curnode->set_beta_(beta_new);
    curnode.beta_ = beta_t(insert_range(curnode.beta_, othnode.beta_, i));
    curnode.support = false;
    assert(curnode.popcount == curnode.rank1(curnode.size()));
}

/*
template <class Vector>
void WaveletTrie::Node::set_beta_(const Vector &bv) {
    assert(bv.size());
    beta_ = beta_t(bv);
    support = false;
}
*/

size_t WaveletTrie::Node::rank0(const size_t i) {
    if (i == size())
        return i - popcount;
    if (!support) {
        sdsl::util::init_support(rank1_, &beta_);
        support = true;
    }
    return i - rank1_(i);
}

size_t WaveletTrie::Node::rank1(const size_t i) {
    if (i == size())
        return popcount;
    if (!support) {
        sdsl::util::init_support(rank1_, &beta_);
        support = true;
    }
    return rank1_(i);
}

template WaveletTrie::WaveletTrie(std::vector<cpp_int>::iterator&, std::vector<cpp_int>::iterator&, size_t);
template WaveletTrie::WaveletTrie(std::vector<cpp_int>&, size_t);
template WaveletTrie::WaveletTrie(std::vector<cpp_int>&&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::set<pos_t>>::iterator&, std::vector<std::set<pos_t>>::iterator&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::set<pos_t>>&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::set<pos_t>>&&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::vector<pos_t>>::iterator&, std::vector<std::vector<pos_t>>::iterator&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::vector<pos_t>>&, size_t);
template WaveletTrie::WaveletTrie(std::vector<std::vector<pos_t>>&&, size_t);

};
