// ldlt/refine.hpp - iterative refinement of a solve.
//
//   x_0 = A \ b (via the factorization);  repeat:  r = b - A x,  x += A \ r
//
// Each step is one residual (with the ORIGINAL A) and one solve with the
// factors.  If the factorization is that of a nearby matrix A + E (static
// pivoting, rounding), refinement converges to the solution of A x = b as long
// as ||A^{-1} E|| < 1, and the normwise backward error drops to O(eps) in a
// few steps.  `a` must therefore be the full symmetric matrix (both triangles),
// so that a.multiply() computes A x.
#pragma once

#include <limits>
#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/numeric.hpp"
#include "ldlt/validate.hpp"


struct RefineOptions {
    Index max_iterations = 10;
    // Stop as soon as the backward error ||r|| / (||A|| ||x|| + ||b||) is
    // below this.  Default: a few units of machine epsilon.
    Real target_backward_error = 4.0 * std::numeric_limits<Real>::epsilon();
    // Stop (stagnation) when a step does not shrink the backward error by at
    // least this factor.  A step that makes it worse is undone.
    Real min_improvement = 0.5;
};

struct RefineResult {
    Index iterations = 0;                   // refinement steps applied
    std::vector<Real> backward_errors;      // [0] = before refinement, one per step
    bool converged = false;                 // reached target_backward_error
    bool stagnated = false;                 // stopped because progress stalled

    Real initial_backward_error() const { return backward_errors.empty() ? 0.0 : backward_errors.front(); }
    Real final_backward_error() const { return backward_errors.empty() ? 0.0 : backward_errors.back(); }
};

// Refine an existing solution x of A x = b in place.
RefineResult refine(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
                    const RefineOptions& opt = {});

// Solve A x = b with the factorization, then refine.
RefineResult solve_refined(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
                           const RefineOptions& opt = {});
std::vector<Real> solve_refined(const CscMatrix& a, const Numeric& num, const std::vector<Real>& b,
                                RefineResult* result = nullptr, const RefineOptions& opt = {});

// check_solve() before and after refinement, for reports and tests.
struct RefinedSolveCheck {
    SolveCheck before;
    SolveCheck after;
    RefineResult refine;
};
RefinedSolveCheck check_solve_refined(const CscMatrix& a, const Numeric& num,
                                      const RefineOptions& opt = {}, bool random_rhs = true,
                                      unsigned seed = 12345u);

