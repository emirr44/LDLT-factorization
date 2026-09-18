#include "csc.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <utility>


namespace {

std::size_t to_size(Index n) { return static_cast<std::size_t>(n < 0 ? 0 : n); }

}

// -- construction --

CscMatrix::CscMatrix(Index nrows, Index ncols)
    : nrows_(nrows), ncols_(ncols), colptr_(to_size(ncols) + 1, 0) {
    if (nrows < 0 || ncols < 0) {
        throw InvalidMatrix("CscMatrix: negative dimensions");
    }
}

CscMatrix::CscMatrix(Index nrows, Index ncols,
                     std::vector<Index> colptr,
                     std::vector<Index> rowind,
                     std::vector<Real> values)
    : nrows_(nrows), ncols_(ncols),
      colptr_(std::move(colptr)), rowind_(std::move(rowind)), values_(std::move(values)) {
    validate();
}

CscMatrix CscMatrix::from_coo(const CooMatrix& coo) {
    coo.validate();
    const Index n = coo.ncols;
    const Index nz = coo.nnz();

    std::vector<Index> colptr(to_size(n) + 1, 0);
    for (Index k = 0; k < nz; ++k) {
        ++colptr[to_size(coo.col[to_size(k)]) + 1];
    }
    for (Index j = 0; j < n; ++j) {
        colptr[to_size(j) + 1] += colptr[to_size(j)];
    }

    std::vector<Index> next(colptr.begin(), colptr.end() - 1);
    std::vector<Index> rowind(to_size(nz));
    std::vector<Real> values(to_size(nz));
    for (Index k = 0; k < nz; ++k) {
        const std::size_t kk = to_size(k);
        const Index p = next[to_size(coo.col[kk])]++;
        rowind[to_size(p)] = coo.row[kk];
        values[to_size(p)] = coo.val[kk];
    }

    CscMatrix a(coo.nrows, n, std::move(colptr), std::move(rowind), std::move(values));
    a.sort_columns();
    return a;
}

CscMatrix CscMatrix::identity(Index n) {
    if (n < 0) throw InvalidMatrix("CscMatrix::identity: negative size");
    std::vector<Index> colptr(to_size(n) + 1);
    std::vector<Index> rowind(to_size(n));
    std::vector<Real> values(to_size(n), 1.0);
    for (Index j = 0; j <= n; ++j) colptr[to_size(j)] = j;
    for (Index j = 0; j < n; ++j) rowind[to_size(j)] = j;
    return CscMatrix(n, n, std::move(colptr), std::move(rowind), std::move(values));
}

CscMatrix CscMatrix::from_dense(Index nrows, Index ncols,
                                const std::vector<Real>& row_major, Real drop_tol) {
    if (nrows < 0 || ncols < 0) throw InvalidMatrix("from_dense: negative dimensions");
    if (row_major.size() != to_size(nrows) * to_size(ncols)) {
        throw InvalidMatrix("from_dense: array size does not match nrows*ncols");
    }
    CooMatrix coo(nrows, ncols);
    for (Index i = 0; i < nrows; ++i) {
        for (Index j = 0; j < ncols; ++j) {
            const Real v = row_major[to_size(i * ncols + j)];
            if (std::abs(v) > drop_tol) coo.add(i, j, v);
        }
    }
    return from_coo(coo);
}

// -- validity --

bool CscMatrix::is_valid(std::string* why) const {
    auto fail = [&](const std::string& msg) {
        if (why) *why = msg;
        return false;
    };
    if (nrows_ < 0 || ncols_ < 0) return fail("negative dimensions");
    if (colptr_.size() != to_size(ncols_) + 1) return fail("colptr must have ncols+1 entries");
    if (colptr_[0] != 0) return fail("colptr[0] must be 0");
    for (Index j = 0; j < ncols_; ++j) {
        if (colptr_[to_size(j)] > colptr_[to_size(j) + 1]) {
            return fail("colptr must be non-decreasing (column " + std::to_string(j) + ")");
        }
    }
    const Index nz = colptr_[to_size(ncols_)];
    if (rowind_.size() != to_size(nz)) return fail("rowind size does not match colptr[ncols]");
    if (values_.size() != to_size(nz)) return fail("values size does not match colptr[ncols]");
    for (Index p = 0; p < nz; ++p) {
        const Index i = rowind_[to_size(p)];
        if (i < 0 || i >= nrows_) {
            return fail("row index " + std::to_string(i) + " out of range at position " +
                        std::to_string(p));
        }
    }
    return true;
}

void CscMatrix::validate() const {
    std::string why;
    if (!is_valid(&why)) throw InvalidMatrix("CscMatrix: " + why);
}

bool CscMatrix::has_sorted_columns() const {
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)] + 1; p < colptr_[to_size(j) + 1]; ++p) {
            if (rowind_[to_size(p) - 1] >= rowind_[to_size(p)]) return false;
        }
    }
    return true;
}

bool CscMatrix::has_duplicates() const {
    std::vector<Index> mark(to_size(nrows_), -1);
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            Index& m = mark[to_size(rowind_[to_size(p)])];
            if (m == j) return true;
            m = j;
        }
    }
    return false;
}

void CscMatrix::sort_columns() {
    std::vector<std::pair<Index, Real>> tmp;
    Index q = 0;   // write position
    for (Index j = 0; j < ncols_; ++j) {
        const Index p0 = colptr_[to_size(j)];
        const Index p1 = colptr_[to_size(j) + 1];
        tmp.clear();
        for (Index p = p0; p < p1; ++p) {
            tmp.emplace_back(rowind_[to_size(p)], values_[to_size(p)]);
        }
        std::sort(tmp.begin(), tmp.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        colptr_[to_size(j)] = q;
        for (std::size_t k = 0; k < tmp.size(); ++k) {
            if (q > colptr_[to_size(j)] && rowind_[to_size(q) - 1] == tmp[k].first) {
                values_[to_size(q) - 1] += tmp[k].second;   // duplicate: sum
            } else {
                rowind_[to_size(q)] = tmp[k].first;
                values_[to_size(q)] = tmp[k].second;
                ++q;
            }
        }
    }
    colptr_[to_size(ncols_)] = q;
    rowind_.resize(to_size(q));
    values_.resize(to_size(q));
}

// -- structure queries --

bool CscMatrix::is_structurally_symmetric() const {
    if (!is_square()) return false;
    CscMatrix a = *this;
    if (!a.has_sorted_columns()) a.sort_columns();
    const CscMatrix t = a.transpose();
    return a.colptr_ == t.colptr_ && a.rowind_ == t.rowind_;
}

bool CscMatrix::is_symmetric(Real tol) const {
    if (!is_square()) return false;
    CscMatrix a = *this;
    if (!a.has_sorted_columns()) a.sort_columns();
    const CscMatrix t = a.transpose();
    if (a.colptr_ != t.colptr_ || a.rowind_ != t.rowind_) return false;
    for (std::size_t p = 0; p < a.values_.size(); ++p) {
        const Real x = a.values_[p];
        const Real y = t.values_[p];
        const Real scale = std::max(std::abs(x), std::abs(y));
        if (std::abs(x - y) > tol * scale) return false;
    }
    return true;
}

std::vector<Real> CscMatrix::diagonal() const {
    const Index n = std::min(nrows_, ncols_);
    std::vector<Real> d(to_size(n), 0.0);
    for (Index j = 0; j < n; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            if (rowind_[to_size(p)] == j) d[to_size(j)] += values_[to_size(p)];
        }
    }
    return d;
}

Real CscMatrix::get(Index i, Index j) const {
    if (i < 0 || i >= nrows_ || j < 0 || j >= ncols_) {
        throw InvalidMatrix("CscMatrix::get: index out of range");
    }
    Real sum = 0.0;
    for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
        if (rowind_[to_size(p)] == i) sum += values_[to_size(p)];
    }
    return sum;
}

CscMatrix::Stats CscMatrix::stats() const {
    Stats s;
    s.nrows = nrows_;
    s.ncols = ncols_;
    s.nnz = nnz();
    s.min_col_nnz = (ncols_ > 0) ? nnz() : 0;
    for (Index j = 0; j < ncols_; ++j) {
        const Index cnt = colptr_[to_size(j) + 1] - colptr_[to_size(j)];
        s.max_col_nnz = std::max(s.max_col_nnz, cnt);
        s.min_col_nnz = std::min(s.min_col_nnz, cnt);
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            const Index i = rowind_[to_size(p)];
            if (i == j) ++s.nnz_diag;
            else if (i > j) ++s.nnz_lower;
            else ++s.nnz_upper;
            s.bandwidth = std::max(s.bandwidth, std::abs(i - j));
        }
    }
    s.avg_col_nnz = (ncols_ > 0) ? static_cast<double>(nnz()) / static_cast<double>(ncols_) : 0.0;
    const double cells = static_cast<double>(nrows_) * static_cast<double>(ncols_);
    s.density = (cells > 0) ? static_cast<double>(nnz()) / cells : 0.0;
    s.structurally_symmetric = is_structurally_symmetric();
    s.numerically_symmetric = s.structurally_symmetric && is_symmetric(1e-12);
    s.norm_1 = norm_1();
    s.norm_inf = norm_inf();
    s.norm_frobenius = norm_frobenius();
    return s;
}

// -- transformations --

CscMatrix CscMatrix::transpose() const {
    // classic counting-sort transpose (cs_transpose)
    const Index nz = nnz();
    std::vector<Index> colptr(to_size(nrows_) + 1, 0);
    for (Index p = 0; p < nz; ++p) ++colptr[to_size(rowind_[to_size(p)]) + 1];
    for (Index i = 0; i < nrows_; ++i) colptr[to_size(i) + 1] += colptr[to_size(i)];

    std::vector<Index> next(colptr.begin(), colptr.end() - 1);
    std::vector<Index> rowind(to_size(nz));
    std::vector<Real> values(to_size(nz));
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            const Index q = next[to_size(rowind_[to_size(p)])]++;
            rowind[to_size(q)] = j;
            values[to_size(q)] = values_[to_size(p)];
        }
    }
    return CscMatrix(ncols_, nrows_, std::move(colptr), std::move(rowind), std::move(values));
}

CscMatrix CscMatrix::filter_triangle(bool keep_upper, bool include_diag) const {
    std::vector<Index> colptr(to_size(ncols_) + 1, 0);
    std::vector<Index> rowind;
    std::vector<Real> values;
    rowind.reserve(rowind_.size());
    values.reserve(values_.size());
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            const Index i = rowind_[to_size(p)];
            bool keep;
            if (i == j) keep = include_diag;
            else keep = keep_upper ? (i < j) : (i > j);
            if (keep) {
                rowind.push_back(i);
                values.push_back(values_[to_size(p)]);
            }
        }
        colptr[to_size(j) + 1] = static_cast<Index>(rowind.size());
    }
    CscMatrix out(nrows_, ncols_, std::move(colptr), std::move(rowind), std::move(values));
    if (!out.has_sorted_columns()) out.sort_columns();
    return out;
}

CscMatrix CscMatrix::upper_triangle(bool include_diag) const {
    return filter_triangle(true, include_diag);
}

CscMatrix CscMatrix::lower_triangle(bool include_diag) const {
    return filter_triangle(false, include_diag);
}

CscMatrix CscMatrix::symmetrize_from_triangle() const {
    if (!is_square()) {
        throw InvalidMatrix("symmetrize_from_triangle: matrix must be square");
    }
    bool has_lower = false, has_upper = false;
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            const Index i = rowind_[to_size(p)];
            if (i > j) has_lower = true;
            if (i < j) has_upper = true;
        }
    }
    if (has_lower && has_upper) {
        throw InvalidMatrix("symmetrize_from_triangle: entries present in both triangles");
    }
    CooMatrix coo(nrows_, ncols_);
    coo.reserve(2 * nnz());
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            coo.add_symmetric(rowind_[to_size(p)], j, values_[to_size(p)]);
        }
    }
    return from_coo(coo);
}

CscMatrix CscMatrix::to_upper_for_factorization() const {
    if (!is_square()) {
        throw InvalidMatrix("to_upper_for_factorization: matrix must be square");
    }
    return upper_triangle(true);
}

// -- arithmetic --

void CscMatrix::multiply(const Real* x, Real* y) const {
    for (Index i = 0; i < nrows_; ++i) y[i] = 0.0;
    for (Index j = 0; j < ncols_; ++j) {
        const Real xj = x[j];
        if (xj == 0.0) continue;
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            y[rowind_[to_size(p)]] += values_[to_size(p)] * xj;
        }
    }
}

std::vector<Real> CscMatrix::multiply(const std::vector<Real>& x) const {
    if (x.size() != to_size(ncols_)) {
        throw InvalidMatrix("multiply: vector length does not match ncols");
    }
    std::vector<Real> y(to_size(nrows_));
    multiply(x.data(), y.data());
    return y;
}

void CscMatrix::multiply_transpose(const Real* x, Real* y) const {
    for (Index j = 0; j < ncols_; ++j) {
        Real sum = 0.0;
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            sum += values_[to_size(p)] * x[rowind_[to_size(p)]];
        }
        y[j] = sum;
    }
}

Real CscMatrix::norm_1() const {
    Real best = 0.0;
    for (Index j = 0; j < ncols_; ++j) {
        Real sum = 0.0;
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            sum += std::abs(values_[to_size(p)]);
        }
        best = std::max(best, sum);
    }
    return best;
}

Real CscMatrix::norm_inf() const {
    std::vector<Real> rowsum(to_size(nrows_), 0.0);
    for (std::size_t p = 0; p < rowind_.size(); ++p) {
        rowsum[to_size(rowind_[p])] += std::abs(values_[p]);
    }
    Real best = 0.0;
    for (Real r : rowsum) best = std::max(best, r);
    return best;
}

Real CscMatrix::norm_frobenius() const {
    Real sum = 0.0;
    for (Real v : values_) sum += v * v;
    return std::sqrt(sum);
}

// -- output --

std::vector<Real> CscMatrix::to_dense() const {
    std::vector<Real> d(to_size(nrows_) * to_size(ncols_), 0.0);
    for (Index j = 0; j < ncols_; ++j) {
        for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1]; ++p) {
            d[to_size(rowind_[to_size(p)] * ncols_ + j)] += values_[to_size(p)];
        }
    }
    return d;
}

void CscMatrix::print(std::ostream& os, Index dense_limit, Index max_triplets) const {
    os << nrows_ << " x " << ncols_ << " sparse matrix, " << nnz() << " stored entries\n";
    if (nrows_ <= dense_limit && ncols_ <= dense_limit) {
        const std::vector<Real> d = to_dense();
        std::ios old_state(nullptr);
        old_state.copyfmt(os);
        os << std::setprecision(4);
        for (Index i = 0; i < nrows_; ++i) {
            for (Index j = 0; j < ncols_; ++j) {
                const Real v = d[to_size(i * ncols_ + j)];
                os << std::setw(10);
                if (v == 0.0) os << ".";
                else os << v;
            }
            os << '\n';
        }
        os.copyfmt(old_state);
    } else {
        Index shown = 0;
        for (Index j = 0; j < ncols_ && shown < max_triplets; ++j) {
            for (Index p = colptr_[to_size(j)]; p < colptr_[to_size(j) + 1] && shown < max_triplets; ++p) {
                os << "  (" << rowind_[to_size(p)] << ", " << j << ")  "
                   << values_[to_size(p)] << '\n';
                ++shown;
            }
        }
        if (shown < nnz()) os << "  ... " << (nnz() - shown) << " more entries\n";
    }
}

std::string CscMatrix::to_string(Index dense_limit, Index max_triplets) const {
    std::ostringstream oss;
    print(oss, dense_limit, max_triplets);
    return oss.str();
}

bool operator==(const CscMatrix& a, const CscMatrix& b) {
    return a.nrows() == b.nrows() && a.ncols() == b.ncols() &&
           a.colptr() == b.colptr() && a.rowind() == b.rowind() && a.values() == b.values();
}

std::ostream& operator<<(std::ostream& os, const CscMatrix& a) {
    a.print(os);
    return os;
}

std::ostream& operator<<(std::ostream& os, const CscMatrix::Stats& s) {
    os << "dimensions            : " << s.nrows << " x " << s.ncols << '\n'
       << "stored entries        : " << s.nnz << '\n'
       << "  diagonal            : " << s.nnz_diag << '\n'
       << "  strictly lower      : " << s.nnz_lower << '\n'
       << "  strictly upper      : " << s.nnz_upper << '\n'
       << "entries per column    : min " << s.min_col_nnz << ", max " << s.max_col_nnz
       << ", avg " << std::setprecision(3) << s.avg_col_nnz << '\n'
       << "bandwidth             : " << s.bandwidth << '\n'
       << "density               : " << std::setprecision(4) << s.density * 100.0 << " %\n"
       << "structurally symmetric: " << (s.structurally_symmetric ? "yes" : "no") << '\n'
       << "numerically symmetric : " << (s.numerically_symmetric ? "yes" : "no") << '\n'
       << std::setprecision(6)
       << "norm_1                : " << s.norm_1 << '\n'
       << "norm_inf              : " << s.norm_inf << '\n'
       << "norm_frobenius        : " << s.norm_frobenius << '\n';
    return os;
}

