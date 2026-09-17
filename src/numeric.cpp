#include "ldlt/numeric.hpp"

#include <algorithm>
#include <cmath>
#include <string>


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }

// The numeric kernel reads A(i,k) with Pinv[i] <= Pinv[k], i.e. the upper
// triangle of P A P'.  When a permutation is present an entry that is upper
// in A may be lower in P A P', and is then only found from the other side of
// the diagonal.  So with a permutation A must store both triangles.  (This is
// an inherited limitation of ldl_numeric; Symbolic::analyze_fast does not
// have it, but Numeric does.)
void check_storage_for_perm(const CscMatrix& a, const Symbolic& s, const char* who) {
    if (!s.has_perm()) return;
    bool has_lower = false, has_upper = false;
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    for (Index j = 0; j < a.ncols() && !(has_lower && has_upper); ++j) {
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i > j) has_lower = true;
            else if (i < j) has_upper = true;
        }
    }
    if (has_upper && !has_lower) {
        throw InvalidMatrix(std::string(who) + ": with a permutation the matrix must store both "
                            "triangles (upper-only storage loses entries that become lower in PAP')");
    }
}
}  // namespace

Numeric Numeric::factorize(const CscMatrix& a, const Symbolic& s, const NumericOptions& opt) {
    if (!a.is_square() || a.nrows() != s.n()) {
        throw InvalidMatrix("Numeric::factorize: matrix does not match the symbolic analysis");
    }
    check_storage_for_perm(a, s, "Numeric::factorize");
    Numeric num;
    num.symbolic_ = s;
    num.li_.assign(sz(s.nnz_L()), 0);
    num.lx_.assign(sz(s.nnz_L()), 0.0);
    num.d_.assign(sz(s.n()), 0.0);
    num.run(a, opt);
    return num;
}

Numeric::Status Numeric::refactorize(const CscMatrix& a, const NumericOptions& opt) {
    if (!a.is_square() || a.nrows() != symbolic_.n()) {
        throw InvalidMatrix("Numeric::refactorize: matrix does not match the symbolic analysis");
    }
    check_storage_for_perm(a, symbolic_, "Numeric::refactorize");
    return run(a, opt);
}

// Port of ldl_numeric.  Comments follow the original.
Numeric::Status Numeric::run(const CscMatrix& a, const NumericOptions& opt) {
    const Index n = symbolic_.n();
    const Index* Ap = a.colptr().data();
    const Index* Ai = a.rowind().data();
    const Real* Ax = a.values().data();
    const Index* Lp = symbolic_.colptr().data();
    const Index* parent = symbolic_.parent().data();
    const Index* P = symbolic_.has_perm() ? symbolic_.perm().data() : nullptr;
    const Index* Pinv = symbolic_.has_perm() ? symbolic_.perm_inv().data() : nullptr;
    Index* Li = li_.data();
    Real* Lx = lx_.data();
    Real* D = d_.data();

    std::vector<Real> y_storage(sz(n), 0.0);
    std::vector<Index> pattern_storage(sz(n)), flag_storage(sz(n)), lnz_storage(sz(n), 0);
    Real* Y = y_storage.data();
    Index* pattern = pattern_storage.data();
    Index* flag = flag_storage.data();
    Index* lnz = lnz_storage.data();

    status_ = Status::Ok;
    failed_pivot_ = -1;
    num_perturbed_ = 0;
    max_perturbation_ = 0.0;
    pivot_threshold_ = (opt.static_pivot_relative > 0.0) ? opt.static_pivot_relative * a.norm_inf() : 0.0;
    const Real threshold = pivot_threshold_;

    for (Index k = 0; k < n; ++k) {
        // compute nonzero Pattern of kth row of L, in topological order
        Y[k] = 0.0;        // Y(0:k) is now all zero
        Index top = n;     // stack for pattern is empty
        flag[k] = k;       // mark node k as visited
        lnz[k] = 0;        // count of nonzeros in column k of L
        const Index kk = P ? P[k] : k;   // kth original, or permuted, column
        const Index p2 = Ap[kk + 1];
        for (Index p = Ap[kk]; p < p2; ++p) {
            Index i = Pinv ? Pinv[Ai[p]] : Ai[p];   // get A(i,k)
            if (i <= k) {
                Y[i] += Ax[p];   // scatter A(i,k) into Y (sum duplicates)
                Index len = 0;
                for (; flag[i] != k; i = parent[i]) {
                    pattern[len++] = i;   // L(k,i) is nonzero
                    flag[i] = k;          // mark i as visited
                }
                while (len > 0) pattern[--top] = pattern[--len];
            }
        }
        // compute numerical values kth row of L (a sparse triangular solve)
        D[k] = Y[k];   // get D(k,k) and clear Y(k)
        Y[k] = 0.0;
        for (; top < n; ++top) {
            const Index i = pattern[top];   // Pattern[top:n-1] is pattern of L(k,:)
            const Real yi = Y[i];           // get and clear Y(i)
            Y[i] = 0.0;
            const Index q2 = Lp[i] + lnz[i];
            for (Index q = Lp[i]; q < q2; ++q) {
                Y[Li[q]] -= Lx[q] * yi;
            }
            const Real l_ki = yi / D[i];   // the nonzero entry L(k,i)
            D[k] -= l_ki * yi;
            Li[q2] = k;                    // store L(k,i) in column form of L
            Lx[q2] = l_ki;
            lnz[i]++;                      // increment count of nonzeros in col i
        }
        if (threshold > 0.0 && std::abs(D[k]) < threshold) {
            // static pivoting: perturb instead of failing (sign kept, + for 0)
            const Real old = D[k];
            D[k] = (old < 0.0) ? -threshold : threshold;
            ++num_perturbed_;
            max_perturbation_ = std::max(max_perturbation_, std::abs(D[k] - old));
        } else if (std::abs(D[k]) <= opt.zero_pivot_tolerance) {   // failure, D(k,k) is zero
            status_ = Status::ZeroPivot;
            failed_pivot_ = k;
            return status_;
        }
    }
    return status_;   // success, diagonal of D is all nonzero
}

CscMatrix Numeric::L() const {
    return CscMatrix(n(), n(), symbolic_.colptr(), li_, lx_);
}

namespace {
Index valid_pivots(Numeric::Status st, Index failed, Index n) {
    return st == Numeric::Status::Ok ? n : failed + 1;
}
}  // namespace

Index Numeric::num_positive_pivots() const {
    const Index m = valid_pivots(status_, failed_pivot_, n());
    Index c = 0;
    for (Index k = 0; k < m; ++k) if (d_[sz(k)] > 0.0) ++c;
    return c;
}

Index Numeric::num_negative_pivots() const {
    const Index m = valid_pivots(status_, failed_pivot_, n());
    Index c = 0;
    for (Index k = 0; k < m; ++k) if (d_[sz(k)] < 0.0) ++c;
    return c;
}

Real Numeric::min_abs_pivot() const {
    const Index m = valid_pivots(status_, failed_pivot_, n());
    if (m == 0) return 0.0;
    Real best = std::abs(d_[0]);
    for (Index k = 1; k < m; ++k) best = std::min(best, std::abs(d_[sz(k)]));
    return best;
}

Real Numeric::max_abs_pivot() const {
    const Index m = valid_pivots(status_, failed_pivot_, n());
    Real best = 0.0;
    for (Index k = 0; k < m; ++k) best = std::max(best, std::abs(d_[sz(k)]));
    return best;
}

// ---------------------------------------------------------------------------
// solves (ports of ldl_lsolve, ldl_dsolve, ldl_ltsolve, ldl_perm, ldl_permt)
// ---------------------------------------------------------------------------

void Numeric::lsolve(Real* x) const {
    const Index n = this->n();
    const Index* Lp = symbolic_.colptr().data();
    for (Index j = 0; j < n; ++j) {
        const Index p2 = Lp[j + 1];
        const Real xj = x[j];
        for (Index p = Lp[j]; p < p2; ++p) {
            x[li_[sz(p)]] -= lx_[sz(p)] * xj;
        }
    }
}

void Numeric::dsolve(Real* x) const {
    const Index n = this->n();
    for (Index j = 0; j < n; ++j) x[j] /= d_[sz(j)];
}

void Numeric::ltsolve(Real* x) const {
    const Index n = this->n();
    const Index* Lp = symbolic_.colptr().data();
    for (Index j = n - 1; j >= 0; --j) {
        const Index p2 = Lp[j + 1];
        Real xj = x[j];
        for (Index p = Lp[j]; p < p2; ++p) {
            xj -= lx_[sz(p)] * x[li_[sz(p)]];
        }
        x[j] = xj;
    }
}

void Numeric::permute(const Real* b, Real* x) const {
    const Index n = this->n();
    const auto& P = symbolic_.perm();
    for (Index j = 0; j < n; ++j) x[j] = b[P[sz(j)]];
}

void Numeric::permute_back(const Real* b, Real* x) const {
    const Index n = this->n();
    const auto& P = symbolic_.perm();
    for (Index j = 0; j < n; ++j) x[P[sz(j)]] = b[j];
}

void Numeric::solve(const Real* b, Real* x) const {
    if (!ok()) throw Error("Numeric::solve: factorization failed (zero pivot)");
    const Index n = this->n();
    if (symbolic_.has_perm()) {
        std::vector<Real> y(sz(n));
        permute(b, y.data());          // y = P b
        lsolve(y.data());
        dsolve(y.data());
        ltsolve(y.data());
        permute_back(y.data(), x);     // x = P' y
    } else {
        if (x != b) std::copy(b, b + n, x);
        lsolve(x);
        dsolve(x);
        ltsolve(x);
    }
}

std::vector<Real> Numeric::solve(const std::vector<Real>& b) const {
    if (static_cast<Index>(b.size()) != n()) {
        throw InvalidMatrix("Numeric::solve: right-hand side has wrong length");
    }
    std::vector<Real> x(b.size());
    solve(b.data(), x.data());
    return x;
}

void Numeric::solve_inplace(std::vector<Real>& x) const {
    if (static_cast<Index>(x.size()) != n()) {
        throw InvalidMatrix("Numeric::solve_inplace: vector has wrong length");
    }
    solve(x.data(), x.data());
}

