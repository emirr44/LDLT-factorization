// symbolic analysis for the LDL' factorization.
//
// Port of ldl_symbolic from Tim Davis's LDL package (LGPL-2.1+), wrapped in a value type that owns its arrays

#pragma once

#include <vector>

#include "csc.hpp"
#include "types.hpp"


class Symbolic {
public:
    Symbolic() = default;

    // Analyze A (must be square)

    // analyze():      port of ldl_symbolic, O(nnz(L)) path following.
    // analyze_fast(): elimination tree + postorder + column counts from CSparse, nearly O(nnz(A)).Produces identical results.
    static Symbolic analyze(const CscMatrix& a, const std::vector<Index>& perm = {});
    static Symbolic analyze_fast(const CscMatrix& a, const std::vector<Index>& perm = {});

    Index n() const { return n_; }
    bool has_perm() const { return !perm_.empty(); }

    const std::vector<Index>& parent() const { return parent_; }
    const std::vector<Index>& col_counts() const { return col_counts_; }
    const std::vector<Index>& colptr() const { return colptr_; }
    const std::vector<Index>& perm() const { return perm_; }          // P
    const std::vector<Index>& perm_inv() const { return perm_inv_; }  // Pinv
    const std::vector<Index>& post() const { return post_; }

    std::vector<Index> postordered_perm() const;

    Index nnz_L() const { return n_ == 0 ? 0 : colptr_.back(); }
    double flops() const;
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

bool is_valid_perm(Index n, const std::vector<Index>& p);

