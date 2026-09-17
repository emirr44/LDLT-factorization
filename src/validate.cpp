#include "ldlt/validate.hpp"

#include <algorithm>
#include <cmath>
#include <random>


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }

Real inf_norm(const Real* v, Index n) {
    Real m = 0.0;
    for (Index i = 0; i < n; ++i) m = std::max(m, std::abs(v[i]));
    return m;
}
}  // namespace

Real Residual::relative() const {
    return b_norm > 0.0 ? r_norm / b_norm : r_norm;
}

Real Residual::backward_error() const {
    const Real denom = a_norm * x_norm + b_norm;
    return denom > 0.0 ? r_norm / denom : r_norm;
}

Residual compute_residual(const CscMatrix& a, const Real* x, const Real* b) {
    if (!a.is_square()) throw InvalidMatrix("compute_residual: matrix must be square");
    const Index n = a.nrows();
    std::vector<Real> ax(sz(n));
    a.multiply(x, ax.data());
    Residual r;
    for (Index i = 0; i < n; ++i) r.r_norm = std::max(r.r_norm, std::abs(b[i] - ax[sz(i)]));
    r.b_norm = inf_norm(b, n);
    r.x_norm = inf_norm(x, n);
    r.a_norm = a.norm_inf();
    return r;
}

Residual compute_residual(const CscMatrix& a, const std::vector<Real>& x,
                          const std::vector<Real>& b) {
    if (static_cast<Index>(x.size()) != a.ncols() || static_cast<Index>(b.size()) != a.nrows()) {
        throw InvalidMatrix("compute_residual: vector length mismatch");
    }
    return compute_residual(a, x.data(), b.data());
}

Real reconstruction_error(const CscMatrix& a, const Numeric& num, Index max_n) {
    const Index n = num.n();
    if (!a.is_square() || a.nrows() != n) {
        throw InvalidMatrix("reconstruction_error: matrix does not match factorization");
    }
    if (!num.ok()) throw Error("reconstruction_error: factorization did not complete");
    if (n > max_n) throw Error("reconstruction_error: n exceeds max_n (dense accumulator)");

    std::vector<Real> m(sz(n) * sz(n), 0.0);
    const auto& Lp = num.L_colptr();
    const auto& Li = num.L_rowind();
    const auto& Lx = num.L_values();
    const auto& D = num.D();

    // M = L D L', one column of L at a time (unit diagonal added explicitly).
    std::vector<Index> rows;
    std::vector<Real> vals;
    for (Index j = 0; j < n; ++j) {
        rows.clear();
        vals.clear();
        rows.push_back(j);
        vals.push_back(1.0);
        for (Index p = Lp[sz(j)]; p < Lp[sz(j) + 1]; ++p) {
            rows.push_back(Li[sz(p)]);
            vals.push_back(Lx[sz(p)]);
        }
        const Real dj = D[sz(j)];
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const Real lr = vals[r] * dj;
            Real* mrow = &m[sz(rows[r]) * sz(n)];
            for (std::size_t c = 0; c < rows.size(); ++c) {
                mrow[sz(rows[c])] += lr * vals[c];
            }
        }
    }

    // M -= P A P'.  Only the upper triangle of A was factorized, so mirror it
    // here; this makes the check independent of whether `a` holds the full
    // matrix or only its upper triangle.
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    const auto& Ax = a.values();
    const bool has_p = num.symbolic().has_perm();
    const auto& Pinv = num.symbolic().perm_inv();
    Real a_frob = 0.0;
    for (Index j = 0; j < n; ++j) {
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i > j) continue;
            const Real v = Ax[sz(p)];
            const Index r = has_p ? Pinv[sz(i)] : i;
            const Index c = has_p ? Pinv[sz(j)] : j;
            m[sz(r) * sz(n) + sz(c)] -= v;
            a_frob += v * v;
            if (i != j) {
                m[sz(c) * sz(n) + sz(r)] -= v;
                a_frob += v * v;
            }
        }
    }
    Real diff = 0.0;
    for (Real v : m) diff += v * v;
    diff = std::sqrt(diff);
    a_frob = std::sqrt(a_frob);
    return a_frob > 0.0 ? diff / a_frob : diff;
}

SolveCheck check_solve(const CscMatrix& a, const Numeric& num, bool random_rhs, unsigned seed) {
    const Index n = num.n();
    if (!a.is_square() || a.nrows() != n) {
        throw InvalidMatrix("check_solve: matrix does not match factorization");
    }
    std::vector<Real> x_true(sz(n), 1.0);
    if (random_rhs) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<Real> dist(-1.0, 1.0);
        for (Real& v : x_true) v = dist(gen);
    }
    const std::vector<Real> b = a.multiply(x_true);
    const std::vector<Real> x = num.solve(b);

    SolveCheck out;
    out.residual = compute_residual(a, x, b);
    Real err = 0.0;
    for (Index i = 0; i < n; ++i) err = std::max(err, std::abs(x[sz(i)] - x_true[sz(i)]));
    const Real xt = inf_norm(x_true.data(), n);
    out.forward_error = xt > 0.0 ? err / xt : err;
    return out;
}

