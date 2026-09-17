#include "ldlt/etree.hpp"

#include <algorithm>


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }
}  // namespace

// ---------------------------------------------------------------------------
// cs_etree (ata = false)
// ---------------------------------------------------------------------------
std::vector<Index> elimination_tree(const CscMatrix& a) {
    if (!a.is_square()) throw InvalidMatrix("elimination_tree: matrix must be square");
    const Index n = a.nrows();
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    std::vector<Index> parent(sz(n), -1);
    std::vector<Index> ancestor(sz(n), -1);   // path-compressed "virtual" parents
    for (Index k = 0; k < n; ++k) {
        for (Index p = Ap[sz(k)]; p < Ap[sz(k) + 1]; ++p) {
            Index i = Ai[sz(p)];
            // A(i,k) with i < k: walk from i towards the root, stop at k
            while (i != -1 && i < k) {
                const Index inext = ancestor[sz(i)];
                ancestor[sz(i)] = k;                       // path compression
                if (inext == -1) parent[sz(i)] = k;        // i was a root: k is its parent
                i = inext;
            }
        }
    }
    return parent;
}

// ---------------------------------------------------------------------------
// cs_post (iterative depth-first search of the forest)
// ---------------------------------------------------------------------------
std::vector<Index> postorder(const std::vector<Index>& parent) {
    const Index n = static_cast<Index>(parent.size());
    std::vector<Index> post(sz(n));
    std::vector<Index> head(sz(n), -1), next(sz(n)), stack(sz(n));
    // youngest-child / next-sibling lists, built so children come out in
    // increasing order
    for (Index j = n - 1; j >= 0; --j) {
        const Index p = parent[sz(j)];
        if (p == -1) continue;
        if (p < 0 || p >= n) throw InvalidMatrix("postorder: parent index out of range");
        next[sz(j)] = head[sz(p)];
        head[sz(p)] = j;
    }
    Index k = 0;
    for (Index j = 0; j < n; ++j) {
        if (parent[sz(j)] != -1) continue;   // only start a DFS at roots
        Index top = 0;
        stack[0] = j;
        while (top >= 0) {
            const Index p = stack[sz(top)];
            const Index i = head[sz(p)];
            if (i == -1) {
                --top;                       // all children done: emit p
                post[sz(k++)] = p;
            } else {
                head[sz(p)] = next[sz(i)];   // remove i from the child list
                stack[sz(++top)] = i;        // descend into i
            }
        }
    }
    if (k != n) throw InvalidMatrix("postorder: parent[] is not a forest");
    return post;
}

// ---------------------------------------------------------------------------
// cs_counts (ata = false)
// ---------------------------------------------------------------------------
namespace {

// Determine whether j is a leaf of the ith row subtree (Gilbert/Ng/Peyton).
// Returns the least common ancestor q of j and the previous leaf, and sets
// jleaf: 0 = j is not a leaf, 1 = first leaf, 2 = subsequent leaf.
Index skeleton_leaf(Index i, Index j, const std::vector<Index>& first,
                    std::vector<Index>& maxfirst, std::vector<Index>& prevleaf,
                    std::vector<Index>& ancestor, int& jleaf) {
    jleaf = 0;
    if (i <= j || first[sz(j)] <= maxfirst[sz(i)]) return -1;
    maxfirst[sz(i)] = first[sz(j)];   // update max first[j] seen so far
    const Index jprev = prevleaf[sz(i)];
    prevleaf[sz(i)] = j;
    jleaf = (jprev == -1) ? 1 : 2;    // j is first or subsequent leaf
    if (jleaf == 1) return i;         // if first leaf, q = root of ith subtree
    Index q = jprev;
    while (q != ancestor[sz(q)]) q = ancestor[sz(q)];
    for (Index s = jprev; s != q;) {  // path compression
        const Index sparent = ancestor[sz(s)];
        ancestor[sz(s)] = q;
        s = sparent;
    }
    return q;   // q = least common ancestor (jprev, j)
}

}  // namespace

std::vector<Index> column_counts(const CscMatrix& a, const std::vector<Index>& parent,
                                 const std::vector<Index>& post) {
    if (!a.is_square()) throw InvalidMatrix("column_counts: matrix must be square");
    const Index n = a.nrows();
    if (static_cast<Index>(parent.size()) != n || static_cast<Index>(post.size()) != n) {
        throw InvalidMatrix("column_counts: parent/post size mismatch");
    }
    // Row-wise access to A: column j of A' is row j of A.
    const CscMatrix at = a.transpose();
    const auto& ATp = at.colptr();
    const auto& ATi = at.rowind();

    std::vector<Index> colcount(sz(n));
    std::vector<Index> ancestor(sz(n)), maxfirst(sz(n), -1), prevleaf(sz(n), -1), first(sz(n), -1);
    std::vector<Index>& delta = colcount;

    // first[j] = smallest postorder index of any descendant of j (incl. j)
    for (Index k = 0; k < n; ++k) {
        Index j = post[sz(k)];
        delta[sz(j)] = (first[sz(j)] == -1) ? 1 : 0;   // delta[j] = 1 if j is a leaf
        for (; j != -1 && first[sz(j)] == -1; j = parent[sz(j)]) first[sz(j)] = k;
    }
    for (Index i = 0; i < n; ++i) ancestor[sz(i)] = i;

    for (Index k = 0; k < n; ++k) {
        const Index j = post[sz(k)];
        if (parent[sz(j)] != -1) delta[sz(parent[sz(j)])]--;   // j is not a root
        for (Index p = ATp[sz(j)]; p < ATp[sz(j) + 1]; ++p) {  // for each A(j,i) in row j
            const Index i = ATi[sz(p)];
            int jleaf = 0;
            const Index q = skeleton_leaf(i, j, first, maxfirst, prevleaf, ancestor, jleaf);
            if (jleaf >= 1) delta[sz(j)]++;
            if (jleaf == 2) delta[sz(q)]--;
        }
        if (parent[sz(j)] != -1) ancestor[sz(j)] = parent[sz(j)];
    }
    // sum up delta's of each child
    for (Index j = 0; j < n; ++j) {
        if (parent[sz(j)] != -1) colcount[sz(parent[sz(j)])] += colcount[sz(j)];
    }
    // cs_counts includes the diagonal; we report strictly-below-diagonal counts
    for (Index j = 0; j < n; ++j) colcount[sz(j)] -= 1;
    return colcount;
}

// ---------------------------------------------------------------------------
// cs_symperm
// ---------------------------------------------------------------------------
CscMatrix permute_symmetric_upper(const CscMatrix& a, const std::vector<Index>& pinv) {
    if (!a.is_square()) throw InvalidMatrix("permute_symmetric_upper: matrix must be square");
    const Index n = a.nrows();
    if (static_cast<Index>(pinv.size()) != n) {
        throw InvalidMatrix("permute_symmetric_upper: pinv has wrong length");
    }
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    const auto& Ax = a.values();

    std::vector<Index> colptr(sz(n) + 1, 0);
    for (Index j = 0; j < n; ++j) {
        const Index j2 = pinv[sz(j)];
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i > j) continue;                  // skip lower triangular part of A
            const Index i2 = pinv[sz(i)];
            colptr[sz(std::max(i2, j2)) + 1]++;   // column count of C
        }
    }
    for (Index j = 0; j < n; ++j) colptr[sz(j) + 1] += colptr[sz(j)];

    std::vector<Index> next(colptr.begin(), colptr.end() - 1);
    std::vector<Index> rowind(sz(colptr[sz(n)]));
    std::vector<Real> values(sz(colptr[sz(n)]));
    for (Index j = 0; j < n; ++j) {
        const Index j2 = pinv[sz(j)];
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i > j) continue;
            const Index i2 = pinv[sz(i)];
            const Index q = next[sz(std::max(i2, j2))]++;
            rowind[sz(q)] = std::min(i2, j2);
            values[sz(q)] = Ax[sz(p)];
        }
    }
    return CscMatrix(n, n, std::move(colptr), std::move(rowind), std::move(values));
}

std::vector<Index> invert_perm(const std::vector<Index>& p) {
    std::vector<Index> pinv(p.size());
    for (Index k = 0; k < static_cast<Index>(p.size()); ++k) pinv[sz(p[sz(k)])] = k;
    return pinv;
}

std::vector<Index> compose_perm(const std::vector<Index>& p, const std::vector<Index>& post) {
    if (p.empty()) return post;
    if (p.size() != post.size()) throw InvalidMatrix("compose_perm: size mismatch");
    std::vector<Index> q(p.size());
    for (std::size_t k = 0; k < p.size(); ++k) q[k] = p[sz(post[k])];
    return q;
}

