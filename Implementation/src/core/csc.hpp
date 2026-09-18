// csc - compressed sparse column matrix

// storage layer used by the symbolic and numeric factorization.

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "coo.hpp"
#include "types.hpp"


class CscMatrix {
public:
    struct Stats {
        Index nrows = 0;
        Index ncols = 0;
        Index nnz = 0;
        Index nnz_diag = 0;        // stored entries with i == j
        Index nnz_lower = 0;       // stored entries with i >  j
        Index nnz_upper = 0;       // stored entries with i <  j
        Index max_col_nnz = 0;
        Index min_col_nnz = 0;
        double avg_col_nnz = 0.0;
        Index bandwidth = 0;       // max |i - j| over stored entries
        double density = 0.0;      // nnz / (nrows * ncols)
        bool structurally_symmetric = false;
        bool numerically_symmetric = false;
        Real norm_1 = 0.0;
        Real norm_inf = 0.0;
        Real norm_frobenius = 0.0;
    };

    CscMatrix() = default;

    CscMatrix(Index nrows, Index ncols);
    CscMatrix(Index nrows, Index ncols,
              std::vector<Index> colptr,
              std::vector<Index> rowind,
              std::vector<Real> values);

    // coo -> csc
    static CscMatrix from_coo(const CooMatrix& coo);
    // n x n identity
    static CscMatrix identity(Index n);

    static CscMatrix from_dense(Index nrows, Index ncols,
                                const std::vector<Real>& row_major,
                                Real drop_tol = 0.0);

    // -- shape --

    Index nrows() const { return nrows_; }
    Index ncols() const { return ncols_; }
    Index nnz() const { return static_cast<Index>(rowind_.size()); }
    bool is_square() const { return nrows_ == ncols_; }
    bool empty() const { return nrows_ == 0 && ncols_ == 0; }

    const std::vector<Index>& colptr() const { return colptr_; }
    const std::vector<Index>& rowind() const { return rowind_; }
    const std::vector<Real>& values() const { return values_; }
    std::vector<Real>& values() { return values_; }  

    // -- validity --

    bool is_valid(std::string* why = nullptr) const;
    void validate() const;
    bool has_sorted_columns() const;
    bool has_duplicates() const;

    void sort_columns();

    // -- structure queries --

    bool is_structurally_symmetric() const;
    // structurally symmetric and |a_ij - a_ji| <= tol * max(|a_ij|, |a_ji|)
    bool is_symmetric(Real tol = 0.0) const;
    std::vector<Real> diagonal() const;
    Real get(Index i, Index j) const;
    Stats stats() const;

    // -- transformations --

    CscMatrix transpose() const;

    CscMatrix upper_triangle(bool include_diag = true) const;
    CscMatrix lower_triangle(bool include_diag = true) const;

    CscMatrix symmetrize_from_triangle() const;
    CscMatrix to_upper_for_factorization() const;

    // -- arithmetic --
    void multiply(const Real* x, Real* y) const;
    std::vector<Real> multiply(const std::vector<Real>& x) const;
    void multiply_transpose(const Real* x, Real* y) const;
    Real norm_1() const;       
    Real norm_inf() const;
    Real norm_frobenius() const;

    // -- output --
    std::vector<Real> to_dense() const;
    void print(std::ostream& os, Index dense_limit = 20,
               Index max_triplets = 50) const;
    std::string to_string(Index dense_limit = 20, Index max_triplets = 50) const;

private:
    Index nrows_ = 0;
    Index ncols_ = 0;
    std::vector<Index> colptr_;
    std::vector<Index> rowind_;
    std::vector<Real> values_;

    CscMatrix filter_triangle(bool keep_upper, bool include_diag) const;
};

bool operator==(const CscMatrix& a, const CscMatrix& b);
inline bool operator!=(const CscMatrix& a, const CscMatrix& b) { return !(a == b); }

std::ostream& operator<<(std::ostream& os, const CscMatrix& a);
std::ostream& operator<<(std::ostream& os, const CscMatrix::Stats& s);

