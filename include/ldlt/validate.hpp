// ldlt/validate.hpp - residual-based validation of a factorization / solve.
#pragma once

#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/numeric.hpp"
#include "ldlt/types.hpp"


// Norms of r = b - A x, all in the infinity norm.
struct Residual {
    Real r_norm = 0.0;     // ||b - A x||
    Real b_norm = 0.0;     // ||b||
    Real x_norm = 0.0;     // ||x||
    Real a_norm = 0.0;     // ||A||_inf

    // ||r|| / ||b||.  Meaningful when b != 0.
    Real relative() const;
    // Normwise backward error  ||r|| / (||A|| ||x|| + ||b||).  For a stable
    // solver this is O(machine epsilon) regardless of conditioning; it is the
    // number to report.
    Real backward_error() const;
};

Residual compute_residual(const CscMatrix& a, const Real* x, const Real* b);
Residual compute_residual(const CscMatrix& a, const std::vector<Real>& x,
                          const std::vector<Real>& b);

// ||P A P' - L D L'||_F / ||A||_F, formed explicitly.  Uses an n-by-n dense
// accumulator, so it is limited to n <= max_n (throws Error above that).
// Intended for tests and small matrices.
Real reconstruction_error(const CscMatrix& a, const Numeric& num, Index max_n = 3000);

// Convenience for the CLI / tests: pick x_true (all ones, or deterministic
// pseudo-random in [-1, 1]), form b = A x_true, solve, and report the
// residual plus the forward error ||x - x_true|| / ||x_true||.
struct SolveCheck {
    Residual residual;
    Real forward_error = 0.0;
};
SolveCheck check_solve(const CscMatrix& a, const Numeric& num, bool random_rhs = true,
                       unsigned seed = 12345u);

