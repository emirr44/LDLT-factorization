#include "ldlt/refine.hpp"

#include <algorithm>
#include <cmath>
#include <random>


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }

Real inf_norm(const std::vector<Real>& v) {
    Real m = 0.0;
    for (Real x : v) m = std::max(m, std::abs(x));
    return m;
}
}  // namespace

RefineResult refine(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
                    const RefineOptions& opt) {
    const Index n = num.n();
    if (!a.is_square() || a.nrows() != n) {
        throw InvalidMatrix("refine: matrix does not match factorization");
    }
    if (!num.ok()) throw Error("refine: factorization did not complete");

    const Real a_norm = a.norm_inf();
    Real b_norm = 0.0;
    for (Index i = 0; i < n; ++i) b_norm = std::max(b_norm, std::abs(b[i]));

    std::vector<Real> r(sz(n)), dx(sz(n)), x_prev(sz(n));
    auto residual = [&]() {
        a.multiply(x, r.data());
        Real r_norm = 0.0, x_norm = 0.0;
        for (Index i = 0; i < n; ++i) {
            r[sz(i)] = b[i] - r[sz(i)];
            r_norm = std::max(r_norm, std::abs(r[sz(i)]));
            x_norm = std::max(x_norm, std::abs(x[i]));
        }
        const Real denom = a_norm * x_norm + b_norm;
        return denom > 0.0 ? r_norm / denom : r_norm;
    };

    RefineResult res;
    Real bwd = residual();
    res.backward_errors.push_back(bwd);
    if (bwd <= opt.target_backward_error) {
        res.converged = true;
        return res;
    }
    for (Index it = 0; it < opt.max_iterations; ++it) {
        std::copy(x, x + n, x_prev.begin());
        num.solve(r.data(), dx.data());                 // dx = (LDL')^{-1} r
        for (Index i = 0; i < n; ++i) x[i] += dx[sz(i)];
        const Real bwd_new = residual();
        if (bwd_new > bwd) {                             // worse: undo and stop
            std::copy(x_prev.begin(), x_prev.end(), x);
            res.stagnated = true;
            break;
        }
        ++res.iterations;
        res.backward_errors.push_back(bwd_new);
        if (bwd_new <= opt.target_backward_error) {
            res.converged = true;
            bwd = bwd_new;
            break;
        }
        if (bwd_new > opt.min_improvement * bwd) {
            res.stagnated = true;
            bwd = bwd_new;
            break;
        }
        bwd = bwd_new;
    }
    return res;
}

RefineResult solve_refined(const CscMatrix& a, const Numeric& num, const Real* b, Real* x,
                           const RefineOptions& opt) {
    num.solve(b, x);
    return refine(a, num, b, x, opt);
}

std::vector<Real> solve_refined(const CscMatrix& a, const Numeric& num, const std::vector<Real>& b,
                                RefineResult* result, const RefineOptions& opt) {
    if (static_cast<Index>(b.size()) != num.n()) {
        throw InvalidMatrix("solve_refined: right-hand side has wrong length");
    }
    std::vector<Real> x(b.size());
    RefineResult r = solve_refined(a, num, b.data(), x.data(), opt);
    if (result) *result = std::move(r);
    return x;
}

RefinedSolveCheck check_solve_refined(const CscMatrix& a, const Numeric& num,
                                      const RefineOptions& opt, bool random_rhs, unsigned seed) {
    const Index n = num.n();
    if (!a.is_square() || a.nrows() != n) {
        throw InvalidMatrix("check_solve_refined: matrix does not match factorization");
    }
    std::vector<Real> x_true(sz(n), 1.0);
    if (random_rhs) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        for (Real& v : x_true) v = dist(gen);
    }
    const std::vector<Real> b = a.multiply(x_true);
    const Real xt_norm = inf_norm(x_true);
    auto forward = [&](const std::vector<Real>& x) {
        Real e = 0.0;
        for (Index i = 0; i < n; ++i) e = std::max(e, std::abs(x[sz(i)] - x_true[sz(i)]));
        return xt_norm > 0.0 ? e / xt_norm : e;
    };

    RefinedSolveCheck out;
    std::vector<Real> x = num.solve(b);
    out.before.residual = compute_residual(a, x, b);
    out.before.forward_error = forward(x);
    out.refine = refine(a, num, b.data(), x.data(), opt);
    out.after.residual = compute_residual(a, x, b);
    out.after.forward_error = forward(x);
    return out;
}

