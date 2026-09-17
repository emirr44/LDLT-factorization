// ldlt/coo.hpp - coordinate (triplet) sparse matrix.
//
// COO is the "assembly" format: entries can be appended in any order and
// duplicates are allowed (they are summed when converting to CSC).  It is the
// natural target for file readers and for building test matrices by hand.
#pragma once

#include <vector>

#include "ldlt/types.hpp"


struct CooMatrix {
    Index nrows = 0;
    Index ncols = 0;
    std::vector<Index> row;   // 0-based row indices
    std::vector<Index> col;   // 0-based column indices
    std::vector<Real> val;    // numerical values

    CooMatrix() = default;
    CooMatrix(Index nrows_, Index ncols_) : nrows(nrows_), ncols(ncols_) {}

    Index nnz() const { return static_cast<Index>(val.size()); }
    bool is_square() const { return nrows == ncols; }

    void reserve(Index capacity);

    // Append the entry A(i, j) += v.  No range checking here; call validate().
    void add(Index i, Index j, Real v);

    // Append (i, j, v) and, when i != j, also (j, i, v).  Convenient for
    // building symmetric matrices from one triangle.
    void add_symmetric(Index i, Index j, Real v);

    // Throws InvalidMatrix if dimensions are negative, the three arrays differ
    // in length, or any index is out of range.
    void validate() const;
};

