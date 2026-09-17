#include "ldlt/symbolic.hpp"

#include "ldlt/etree.hpp"


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }
}  // namespace

bool is_valid_perm(Index n, const std::vector<Index>& p) {
    if (p.empty()) return true;
    if (static_cast<Index>(p.size()) != n) return false;
    std::vector<char> seen(sz(n), 0);
    for (Index j : p) {
        if (j < 0 || j >= n) return false;
        char& s = seen[sz(j)];
        if (s) return false;
        s = 1;
    }
    return true;
}

namespace {
// Counts stored entries strictly below / above the diagonal.
void count_triangles(const CscMatrix& a, Index& lower, Index& upper) {
    lower = upper = 0;
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    for (Index j = 0; j < a.ncols(); ++j) {
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i > j) ++lower;
            else if (i < j) ++upper;
        }
    }
}
}  // namespace

// Shared validation and permutation setup.
Symbolic Symbolic::prepare(const CscMatrix& a, const std::vector<Index>& perm) {
    if (!a.is_square()) throw InvalidMatrix("Symbolic::analyze: matrix must be square");
    const Index n = a.nrows();
    if (!is_valid_perm(n, perm)) {
        throw InvalidMatrix("Symbolic::analyze: perm is not a permutation of 0..n-1");
    }
    // The factorization reads the upper triangle.  A matrix that stores only
    // its lower triangle would silently factor as a diagonal matrix (this is
    // what happens with an unexpanded symmetric Matrix Market file), so
    // refuse it here.
    Index lower = 0, upper = 0;
    count_triangles(a, lower, upper);
    if (lower > 0 && upper == 0) {
        throw InvalidMatrix("Symbolic::analyze: matrix stores only its lower triangle; "
                            "LDL' reads the upper triangle - pass the full symmetric matrix "
                            "(e.g. CscMatrix::symmetrize_from_triangle())");
    }
    Symbolic s;
    s.n_ = n;
    s.parent_.assign(sz(n), -1);
    s.col_counts_.assign(sz(n), 0);
    s.colptr_.assign(sz(n) + 1, 0);
    if (!perm.empty()) {
        s.perm_ = perm;
        s.perm_inv_ = invert_perm(perm);
    }
    return s;
}

void Symbolic::finish_colptr() {
    for (Index k = 0; k < n_; ++k) {
        colptr_[sz(k) + 1] = colptr_[sz(k)] + col_counts_[sz(k)];
    }
}

// ---------------------------------------------------------------------------
// analyze: port of ldl_symbolic
// ---------------------------------------------------------------------------
Symbolic Symbolic::analyze(const CscMatrix& a, const std::vector<Index>& perm) {
    Symbolic s = prepare(a, perm);
    const Index n = s.n_;
    const bool has_p = !perm.empty();

    // Raw pointers keep the inner loop identical to ldl_symbolic.
    const Index* Ap = a.colptr().data();
    const Index* Ai = a.rowind().data();
    const Index* P = has_p ? s.perm_.data() : nullptr;
    const Index* Pinv = has_p ? s.perm_inv_.data() : nullptr;
    Index* parent = s.parent_.data();
    Index* lnz = s.col_counts_.data();
    std::vector<Index> flag_storage(sz(n));
    Index* flag = flag_storage.data();

    for (Index k = 0; k < n; ++k) {
        // L(k,:) pattern: all nodes reachable in the etree from nz in A(0:k-1,k)
        parent[k] = -1;   // parent of k is not yet known
        flag[k] = k;      // mark node k as visited
        lnz[k] = 0;       // count of nonzeros in column k of L
        const Index kk = P ? P[k] : k;   // kth original, or permuted, column
        const Index p2 = Ap[kk + 1];
        for (Index p = Ap[kk]; p < p2; ++p) {
            // A(i,k) is nonzero (original or permuted A)
            Index i = Pinv ? Pinv[Ai[p]] : Ai[p];
            if (i < k) {
                // follow path from i to root of etree, stop at flagged node
                for (; flag[i] != k; i = parent[i]) {
                    if (parent[i] == -1) parent[i] = k;   // find parent of i if not yet determined
                    lnz[i]++;                             // L(k,i) is nonzero
                    flag[i] = k;                          // mark i as visited
                }
            }
        }
    }
    s.finish_colptr();
    s.post_ = postorder(s.parent_);
    return s;
}

// ---------------------------------------------------------------------------
// analyze_fast: etree + postorder + column counts
// ---------------------------------------------------------------------------
Symbolic Symbolic::analyze_fast(const CscMatrix& a, const std::vector<Index>& perm) {
    Symbolic s = prepare(a, perm);
    // Work on the upper triangle of P A P' (only formed when there is a P).
    CscMatrix permuted;
    const CscMatrix* c = &a;
    if (s.has_perm()) {
        permuted = permute_symmetric_upper(a, s.perm_inv_);
        c = &permuted;
    }
    s.parent_ = elimination_tree(*c);
    s.post_ = postorder(s.parent_);
    s.col_counts_ = column_counts(*c, s.parent_, s.post_);
    s.finish_colptr();
    return s;
}

std::vector<Index> Symbolic::postordered_perm() const {
    return compose_perm(perm_, post_);
}

double Symbolic::flops() const {
    double f = 0.0;
    for (Index c : col_counts_) {
        const double cc = static_cast<double>(c);
        f += cc * (cc + 2.0);
    }
    return f;
}

double Symbolic::fill_ratio(const CscMatrix& a) const {
    // count entries of A on or above the diagonal (the part that is read)
    Index nnz_upper = 0;
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    for (Index j = 0; j < a.ncols(); ++j) {
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            if (Ai[sz(p)] <= j) ++nnz_upper;
        }
    }
    if (nnz_upper == 0) return 0.0;
    return static_cast<double>(nnz_L() + n_) / static_cast<double>(nnz_upper);
}

