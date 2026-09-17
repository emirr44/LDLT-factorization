// ldlt/csc.hpp - compressed sparse column matrix.
//
// This is the storage layer used by the symbolic and numeric factorization.
// Column j occupies positions colptr[j] .. colptr[j+1]-1 of rowind/values.
// Invariants maintained by every constructor and operation in this class:
//   * colptr has size ncols+1, colptr[0] == 0, and is non-decreasing
//   * rowind and values have size colptr[ncols]
//   * every row index is in [0, nrows)
// Sorted, duplicate-free columns are additionally guaranteed after from_coo(),
// transpose(), sort_columns(), and every triangle/symmetrize operation.
// The factorization code will rely on that, so prefer those entry points.
#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "ldlt/coo.hpp"
#include "ldlt/types.hpp"


class CscMatrix {
public:
    // Summary statistics, mainly for the command line "info" output.
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
        bool numerically_symmetric = false;   // within a small relative tolerance
        Real norm_1 = 0.0;
        Real norm_inf = 0.0;
        Real norm_frobenius = 0.0;
    };

    CscMatrix() = default;

    // Empty nrows x ncols matrix (no stored entries).
    CscMatrix(Index nrows, Index ncols);

    // Takes ownership of the three arrays.  Validates the invariants above and
    // throws InvalidMatrix if they fail.  Columns are not sorted automatically.
    CscMatrix(Index nrows, Index ncols,
              std::vector<Index> colptr,
              std::vector<Index> rowind,
              std::vector<Real> values);

    // Conversion from triplet form.  Duplicate (i, j) pairs are summed and the
    // rows within each column come out sorted ascending.
    static CscMatrix from_coo(const CooMatrix& coo);

    // n x n identity.
    static CscMatrix identity(Index n);

    // Build directly from a dense row-major array (tests, small examples).
    // Entries with |a_ij| <= drop_tol are not stored.
    static CscMatrix from_dense(Index nrows, Index ncols,
                                const std::vector<Real>& row_major,
                                Real drop_tol = 0.0);

    // ---- shape and raw access -------------------------------------------
    Index nrows() const { return nrows_; }
    Index ncols() const { return ncols_; }
    Index nnz() const { return static_cast<Index>(rowind_.size()); }
    bool is_square() const { return nrows_ == ncols_; }
    bool empty() const { return nrows_ == 0 && ncols_ == 0; }

    const std::vector<Index>& colptr() const { return colptr_; }
    const std::vector<Index>& rowind() const { return rowind_; }
    const std::vector<Real>& values() const { return values_; }
    std::vector<Real>& values() { return values_; }   // numeric refactorization

    // ---- validity --------------------------------------------------------
    // Returns true if the structural invariants hold.  When `why` is non-null
    // it receives a human-readable reason on failure.
    bool is_valid(std::string* why = nullptr) const;
    void validate() const;   // throws InvalidMatrix
    bool has_sorted_columns() const;   // strictly increasing rows, so no dups
    bool has_duplicates() const;

    // Sort each column by row index and sum duplicate entries.
    void sort_columns();

    // ---- structure queries --------------------------------------------------
    bool is_structurally_symmetric() const;
    // Structurally symmetric and |a_ij - a_ji| <= tol * max(|a_ij|, |a_ji|).
    bool is_symmetric(Real tol = 0.0) const;
    // Diagonal as a dense vector of length min(nrows, ncols); missing => 0.
    std::vector<Real> diagonal() const;
    // A(i, j) by linear search of column j.  O(nnz in column j).
    Real get(Index i, Index j) const;
    Stats stats() const;

    // ---- transformations ----------------------------------------------------
    CscMatrix transpose() const;
    // Keep only entries with i <= j (or i < j when include_diag is false).
    CscMatrix upper_triangle(bool include_diag = true) const;
    // Keep only entries with i >= j (or i > j when include_diag is false).
    CscMatrix lower_triangle(bool include_diag = true) const;
    // Given a matrix that stores one triangle of a symmetric matrix (either
    // one), return the full symmetric matrix T + T' - diag(T).  Throws if
    // entries exist strictly on both sides of the diagonal.
    CscMatrix symmetrize_from_triangle() const;
    // Given a full symmetric matrix return only the upper triangle, which is
    // the part LDL' reads.  Equivalent to upper_triangle() but checks that the
    // matrix is square.
    CscMatrix to_upper_for_factorization() const;

    // ---- arithmetic -------------------------------------------------------------
    // y = A x   (x has ncols entries, y has nrows entries)
    void multiply(const Real* x, Real* y) const;
    std::vector<Real> multiply(const std::vector<Real>& x) const;
    // y = A' x  (x has nrows entries, y has ncols entries)
    void multiply_transpose(const Real* x, Real* y) const;
    Real norm_1() const;          // max column sum of |a_ij|
    Real norm_inf() const;        // max row sum of |a_ij|
    Real norm_frobenius() const;

    // ---- output -----------------------------------------------------------------
    // Dense copy, row-major, nrows*ncols entries.  For small matrices only.
    std::vector<Real> to_dense() const;
    // Pretty-print.  Dense grid when both dimensions are <= dense_limit,
    // otherwise a triplet listing (at most max_triplets entries).
    void print(std::ostream& os, Index dense_limit = 20,
               Index max_triplets = 50) const;
    std::string to_string(Index dense_limit = 20, Index max_triplets = 50) const;

private:
    Index nrows_ = 0;
    Index ncols_ = 0;
    std::vector<Index> colptr_;
    std::vector<Index> rowind_;
    std::vector<Real> values_;

    // keep_upper: entries with i <= j (strict when !include_diag)
    CscMatrix filter_triangle(bool keep_upper, bool include_diag) const;
};

// Equality of shape, structure and values (exact).
bool operator==(const CscMatrix& a, const CscMatrix& b);
inline bool operator!=(const CscMatrix& a, const CscMatrix& b) { return !(a == b); }

std::ostream& operator<<(std::ostream& os, const CscMatrix& a);
std::ostream& operator<<(std::ostream& os, const CscMatrix::Stats& s);

