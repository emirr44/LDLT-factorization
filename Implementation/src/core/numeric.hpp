// numeric LDL' factorization and triangular solves.

#pragma once

#include <vector>

#include "csc.hpp"
#include "symbolic.hpp"
#include "types.hpp"


struct NumericOptions {
    /*
    a pivot d with |d| <= zero_pivot_tolerance is treated as zero and the
    factorization stops (status ZeroPivot). The original ldl.c uses an
    exact test, i.e. 0.0
    */
    Real zero_pivot_tolerance = 0.0;

    /*
    static pivoting.  If > 0, every pivot with |d| < threshold, where
    threshold = static_pivot_relative * ||A||_inf, is replaced by
    +-threshold (sign kept, + for an exact zero). The factorization then
    never fails, but it is the exact factorization of a slightly perturbed
    matrix A + E with ||E|| <= threshold per perturbed pivot, so the solve
    must be followed by iterative refinement to recover
    a small backward error for A itself
    */
    Real static_pivot_relative = 0.0;
};


class Numeric {
public:
    enum class Status { Ok, ZeroPivot };

    Numeric() = default;

    static Numeric factorize(const CscMatrix& a, const Symbolic& s,
                             const NumericOptions& opt = {});

    // Recompute L and D for a matrix with the SAME pattern (and permutation) as the one this object was created from
    Status refactorize(const CscMatrix& a, const NumericOptions& opt = {});

    Status status() const { return status_; }
    bool ok() const { return status_ == Status::Ok; }
    Index failed_pivot() const { return failed_pivot_; }

    Index n() const { return symbolic_.n(); }
    const Symbolic& symbolic() const { return symbolic_; }

    CscMatrix L() const;
    const std::vector<Index>& L_colptr() const { return symbolic_.colptr(); }
    const std::vector<Index>& L_rowind() const { return li_; }
    const std::vector<Real>& L_values() const { return lx_; }
    const std::vector<Real>& D() const { return d_; }

    Index nnz_L() const { return symbolic_.nnz_L(); }
    Index num_positive_pivots() const;
    Index num_negative_pivots() const;
    Real min_abs_pivot() const;
    Real max_abs_pivot() const;

    Index num_perturbed_pivots() const { return num_perturbed_; }
    Real max_perturbation() const { return max_perturbation_; }
    Real pivot_threshold() const { return pivot_threshold_; }

    // solve A x = b for one right-hand side
    void solve(const Real* b, Real* x) const;
    std::vector<Real> solve(const std::vector<Real>& b) const;
    void solve_inplace(std::vector<Real>& x) const;   // x is b on input

    void lsolve(Real* x) const;    // x := L  \ x
    void dsolve(Real* x) const;    // x := D  \ x
    void ltsolve(Real* x) const;   // x := L' \ x
    void permute(const Real* b, Real* x) const;         // x = P b   (x[k] = b[P[k]])
    void permute_back(const Real* b, Real* x) const;    // x = P' b  (x[P[k]] = b[k])

private:
    Symbolic symbolic_;
    std::vector<Index> li_;   // row indices of L
    std::vector<Real> lx_;    // values of L
    std::vector<Real> d_;     // diagonal of D
    Status status_ = Status::Ok;
    Index failed_pivot_ = -1;
    Index num_perturbed_ = 0;
    Real max_perturbation_ = 0.0;
    Real pivot_threshold_ = 0.0;

    Status run(const CscMatrix& a, const NumericOptions& opt);
};

