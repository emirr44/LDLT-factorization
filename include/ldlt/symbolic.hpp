// ldlt/symbolic.hpp - symbolic analysis for the LDL' factorization.
//
// Port of ldl_symbolic from Tim Davis's LDL package (LGPL-2.1+), wrapped in a
// value type that owns its arrays.  Given the pattern of the upper triangle of
// a symmetric n-by-n matrix A (and an optional fill-reducing permutation P),
// it computes
//   * the elimination tree        parent[j]  (-1 for a root)
//   * the column counts of L      col_counts[j] = nnz(L(:,j)) below the diagonal
//   * the column pointers of L    colptr[0..n]
//   * P and its inverse           perm[k] = j means row/col j of A is pivot k
//
// Only entries A(i,j) with i < j (in the permuted matrix) are read.  Anything
// in the strictly lower triangle is ignored, so callers may pass the full
// symmetric matrix or only its upper triangle.  Passing only the LOWER
// triangle silently yields a diagonal factorization - see CscMatrix docs.
//
// Cost: O(nnz(L)), because it follows the etree path for every entry of A.
// (CSparse's cs_etree/cs_counts do this in near O(nnz(A)); that is a later
// phase.)
#pragma once

#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/types.hpp"


class Symbolic {
public:
    Symbolic() = default;

    // Analyze A (must be square).  `perm` is either empty (natural ordering)
    // or a permutation of 0..n-1; it is validated and InvalidMatrix is thrown
    // if it is not a permutation.
    //
    // analyze():      port of ldl_symbolic, O(nnz(L)) path following.
    // analyze_fast(): elimination tree + postorder + column counts from
    //                 CSparse, nearly O(nnz(A)).  Produces identical results.
    static Symbolic analyze(const CscMatrix& a, const std::vector<Index>& perm = {});
    static Symbolic analyze_fast(const CscMatrix& a, const std::vector<Index>& perm = {});

    Index n() const { return n_; }
    bool has_perm() const { return !perm_.empty(); }

    const std::vector<Index>& parent() const { return parent_; }
    const std::vector<Index>& col_counts() const { return col_counts_; }
    const std::vector<Index>& colptr() const { return colptr_; }
    const std::vector<Index>& perm() const { return perm_; }          // P
    const std::vector<Index>& perm_inv() const { return perm_inv_; }  // Pinv
    // Postorder of the elimination tree (in the permuted numbering).
    const std::vector<Index>& post() const { return post_; }
    // The ordering P composed with the postorder: same fill, but the etree
    // becomes postordered, which improves locality and is what supernodal
    // methods require.  Returns post() itself when there is no permutation.
    std::vector<Index> postordered_perm() const;

    // nnz of L strictly below the diagonal (= colptr[n]).
    Index nnz_L() const { return n_ == 0 ? 0 : colptr_.back(); }
    // Floating-point operations of the numeric phase: sum over k of
    // col_counts[k] * (col_counts[k] + 2).
    double flops() const;
    // nnz(L + D) / nnz(upper triangle of A incl. diagonal); needs A again.
    double fill_ratio(const CscMatrix& a) const;

private:
    Index n_ = 0;
    std::vector<Index> parent_;
    std::vector<Index> col_counts_;
    std::vector<Index> colptr_;
    std::vector<Index> post_;
    std::vector<Index> perm_;
    std::vector<Index> perm_inv_;

    static Symbolic prepare(const CscMatrix& a, const std::vector<Index>& perm);
    void finish_colptr();
};

// Returns true if p is a permutation of 0..n-1.  Empty p counts as identity.
bool is_valid_perm(Index n, const std::vector<Index>& p);

