// iterative refinement of a solve.

#pragma once

#include <limits>
#include <vector>

#include "csc.hpp"
#include "numeric.hpp"
#include "validate.hpp"


struct RefineOptions {
    Index max_iterations = 10;
    // stop as soon as the backward error ||r|| / (||A|| ||x|| + ||b||) is below this
    Real target_backward_error = 4.0 * std::numeric_limits<Real>::epsilon();
    // stop (stagnation) when a step does not shrink the backward error by at least this factor
    Real min_improvement = 0.5;
};

struct RefineResult {
    Index iterations = 0;                   
    std::vector<Real> backward_errors;
    bool converged = false;              
    bool stagnated = false;

    Real initial_backward_error() const { return backward_errors.empty() ? 0.0 : backward_errors.front(); }
    Real final_backward_error() const { return backward_errors.empty() ? 0.0 : backward_errors.back(); }
};

RefineResult refine(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
        const RefineOptions& opt = {});

RefineResult solve_refined(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
        const RefineOptions& opt = {});
std::vector<Real> solve_refined(const CscMatrix& a, const Numeric& num, const std::vector<Real>& b,
        RefineResult* result = nullptr, const RefineOptions& opt = {});

struct RefinedSolveCheck {
    SolveCheck before;
    SolveCheck after;
    RefineResult refine;
};
RefinedSolveCheck check_solve_refined(const CscMatrix& a, const Numeric& num,
        const RefineOptions& opt = {}, bool random_rhs = true, unsigned seed = 12345u);

