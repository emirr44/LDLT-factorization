// residual-based validation of a factorization / solve
#pragma once

#include <vector>

#include "csc.hpp"
#include "numeric.hpp"
#include "types.hpp"


struct Residual {
    Real r_norm = 0.0;     // ||b - A x||
    Real b_norm = 0.0;     // ||b||
    Real x_norm = 0.0;     // ||x||
    Real a_norm = 0.0;     // ||A||_inf

    // ||r|| / ||b||.  Meaningful when b != 0.
    Real relative() const;
    // Normwise backward error  ||r|| / (||A|| ||x|| + ||b||)
    Real backward_error() const;
};

Residual compute_residual(const CscMatrix& a, const Real* x, const Real* b);
Residual compute_residual(const CscMatrix& a, const std::vector<Real>& x,
        const std::vector<Real>& b);

// ||P A P' - L D L'||_F / ||A||_F, formed explicitly
Real reconstruction_error(const CscMatrix& a, const Numeric& num, Index max_n = 3000);

struct SolveCheck {
    Residual residual;
    Real forward_error = 0.0;
};
SolveCheck check_solve(const CscMatrix& a, const Numeric& num, bool random_rhs = true,
        unsigned seed = 12345u);

