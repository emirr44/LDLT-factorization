// ldlt/numeric.hpp - numeric LDL' factorization and triangular solves.
//
// Port of ldl_numeric / ldl_lsolve / ldl_dsolve / ldl_ltsolve / ldl_perm /
// ldl_permt from Tim Davis's LDL package (LGPL-2.1+), wrapped in a class.
//
//   P A P' = L D L'
//
// L is unit lower triangular (the unit diagonal is NOT stored), D is diagonal.
// The algorithm is up-looking: row k of L is obtained from a sparse triangular
// solve with L(0:k-1, 0:k-1), following the elimination tree.  No pivoting is
// performed, so the factorization is stable for symmetric positive definite
// matrices and may be inaccurate for indefinite ones (later phase).
#pragma once

#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/symbolic.hpp"
#include "ldlt/types.hpp"


struct NumericOptions {
    // A pivot d with |d| <= zero_pivot_tolerance is treated as zero and the
    // factorization stops (status ZeroPivot).  The original ldl.c uses an
    // exact test, i.e. 0.0.  A relative threshold such as 1e-14 * norm(A) is
    // usually more sensible; the caller decides.
    Real zero_pivot_tolerance = 0.0;

    // Static pivoting.  If > 0, every pivot with |d| < threshold, where
    // threshold = static_pivot_relative * ||A||_inf, is replaced by
    // +-threshold (sign kept, + for an exact zero).  The factorization then
    // never fails, but it is the exact factorization of a slightly perturbed
    // matrix A + E with ||E|| <= threshold per perturbed pivot, so the solve
    // must be followed by iterative refinement (ldlt/refine.hpp) to recover
    // a small backward error for A itself.  sqrt(eps) ~ 1e-8 is the classic
    // choice (Li & Demmel; PARDISO uses the same idea).  0 disables it.
    Real static_pivot_relative = 0.0;
};

class Numeric {
public:
    enum class Status { Ok, ZeroPivot };

    Numeric() = default;

    // Factorize A using a previous symbolic analysis of the same pattern.
    // A must be square with n == s.n().  Only the upper triangle of P A P' is
    // read; therefore, when the analysis carries a permutation, A must store
    // BOTH triangles (throws InvalidMatrix for upper-only storage).  Without
    // a permutation the full matrix or its upper triangle may be passed.
    // Never throws for numerical trouble: check status() afterwards.
    static Numeric factorize(const CscMatrix& a, const Symbolic& s,
                             const NumericOptions& opt = {});

    // Recompute L and D for a matrix with the SAME pattern (and permutation)
    // as the one this object was created from.  Reuses all storage.
    Status refactorize(const CscMatrix& a, const NumericOptions& opt = {});

    Status status() const { return status_; }
    bool ok() const { return status_ == Status::Ok; }
    // Index k of the pivot that was (numerically) zero, or -1 if ok().
    Index failed_pivot() const { return failed_pivot_; }

    Index n() const { return symbolic_.n(); }
    const Symbolic& symbolic() const { return symbolic_; }

    // L as a CSC matrix without its unit diagonal (rows sorted ascending).
    // Only meaningful when ok(); on failure columns >= failed_pivot() are
    // incomplete.
    CscMatrix L() const;
    const std::vector<Index>& L_colptr() const { return symbolic_.colptr(); }
    const std::vector<Index>& L_rowind() const { return li_; }
    const std::vector<Real>& L_values() const { return lx_; }
    const std::vector<Real>& D() const { return d_; }

    // nnz of L below the diagonal (the whole factor when ok()).
    Index nnz_L() const { return symbolic_.nnz_L(); }
    // Inertia information (Sylvester): pivots of each sign.
    Index num_positive_pivots() const;
    Index num_negative_pivots() const;
    Real min_abs_pivot() const;
    Real max_abs_pivot() const;
    // Static pivoting statistics of the last factorization.
    Index num_perturbed_pivots() const { return num_perturbed_; }
    Real max_perturbation() const { return max_perturbation_; }
    Real pivot_threshold() const { return pivot_threshold_; }   // 0 if disabled

    // Solve A x = b for one right-hand side.  Requires ok().
    // Handles the permutation internally:  x = P' (L' \ (D \ (L \ (P b)))).
    void solve(const Real* b, Real* x) const;
    std::vector<Real> solve(const std::vector<Real>& b) const;
    void solve_inplace(std::vector<Real>& x) const;   // x is b on input

    // Building blocks, operating on a vector in the PERMUTED space.
    void lsolve(Real* x) const;    // x := L  \ x
    void dsolve(Real* x) const;    // x := D  \ x
    void ltsolve(Real* x) const;   // x := L' \ x
    void permute(const Real* b, Real* x) const;         // x = P b   (x[k] = b[P[k]])
    void permute_back(const Real* b, Real* x) const;    // x = P' b  (x[P[k]] = b[k])

private:
    Symbolic symbolic_;
    std::vector<Index> li_;   // row indices of L, size colptr[n]
    std::vector<Real> lx_;    // values of L
    std::vector<Real> d_;     // diagonal of D
    Status status_ = Status::Ok;
    Index failed_pivot_ = -1;
    Index num_perturbed_ = 0;
    Real max_perturbation_ = 0.0;
    Real pivot_threshold_ = 0.0;

    Status run(const CscMatrix& a, const NumericOptions& opt);
};

