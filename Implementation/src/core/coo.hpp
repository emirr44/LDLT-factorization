// coordinate (triplet) sparse matrix.

#pragma once

#include <vector>

#include "types.hpp"


struct CooMatrix {
    Index nrows = 0;
    Index ncols = 0;
    std::vector<Index> row;
    std::vector<Index> col;
    std::vector<Real> val;

    CooMatrix() = default;
    CooMatrix(Index nrows_, Index ncols_) : nrows(nrows_), ncols(ncols_) {}

    Index nnz() const { return static_cast<Index>(val.size()); }
    bool is_square() const { return nrows == ncols; }

    void reserve(Index capacity);

    void add(Index i, Index j, Real v);
    void add_symmetric(Index i, Index j, Real v);

    void validate() const;
};

