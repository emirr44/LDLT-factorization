#include "ldlt/coo.hpp"

#include <string>


void CooMatrix::reserve(Index capacity) {
    if (capacity < 0) capacity = 0;
    row.reserve(static_cast<std::size_t>(capacity));
    col.reserve(static_cast<std::size_t>(capacity));
    val.reserve(static_cast<std::size_t>(capacity));
}

void CooMatrix::add(Index i, Index j, Real v) {
    row.push_back(i);
    col.push_back(j);
    val.push_back(v);
}

void CooMatrix::add_symmetric(Index i, Index j, Real v) {
    add(i, j, v);
    if (i != j) add(j, i, v);
}

void CooMatrix::validate() const {
    if (nrows < 0 || ncols < 0) {
        throw InvalidMatrix("CooMatrix: negative dimensions");
    }
    if (row.size() != col.size() || row.size() != val.size()) {
        throw InvalidMatrix("CooMatrix: row/col/val arrays differ in length");
    }
    for (std::size_t k = 0; k < row.size(); ++k) {
        if (row[k] < 0 || row[k] >= nrows) {
            throw InvalidMatrix("CooMatrix: row index " + std::to_string(row[k]) +
                                " out of range at entry " + std::to_string(k));
        }
        if (col[k] < 0 || col[k] >= ncols) {
            throw InvalidMatrix("CooMatrix: column index " + std::to_string(col[k]) +
                                " out of range at entry " + std::to_string(k));
        }
    }
}

