// ldlt/matrix_market.hpp - Matrix Market (.mtx) reader and writer.
//
// Supported: "matrix coordinate {real|integer|pattern} {general|symmetric|
// skew-symmetric}".  Complex and dense ("array") files are rejected with a
// clear IoError.  Indices in the file are 1-based; in memory they are 0-based.
//
// Important: symmetric Matrix Market files store only the LOWER triangle
// (i >= j).  The LDL' factorization reads the UPPER triangle, so a file read
// without expansion and fed straight to the factorization would factor a
// diagonal matrix.  By default read_matrix_market() therefore expands
// symmetric/skew-symmetric files to the full matrix.
#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "ldlt/coo.hpp"
#include "ldlt/csc.hpp"


struct MatrixMarketInfo {
    enum class Field { Double, Integer, Pattern };   // "real" in the file = Double
    enum class Symmetry { General, Symmetric, SkewSymmetric };

    Field field = Field::Double;
    Symmetry symmetry = Symmetry::General;
    Index nrows = 0;
    Index ncols = 0;
    Index nnz_in_file = 0;          // entries listed in the file
    std::vector<std::string> comments;   // '%' lines after the header

    bool is_symmetric_storage() const { return symmetry != Symmetry::General; }
};

// Read into triplet form.  When expand_symmetric is true, symmetric and
// skew-symmetric files produce both (i, j) and (j, i) entries.
CooMatrix read_matrix_market_coo(std::istream& in,
                                 bool expand_symmetric = true,
                                 MatrixMarketInfo* info = nullptr);
CooMatrix read_matrix_market_coo(const std::string& path,
                                 bool expand_symmetric = true,
                                 MatrixMarketInfo* info = nullptr);

// Read directly into CSC (sorted columns, duplicates summed).
CscMatrix read_matrix_market(std::istream& in,
                             bool expand_symmetric = true,
                             MatrixMarketInfo* info = nullptr);
CscMatrix read_matrix_market(const std::string& path,
                             bool expand_symmetric = true,
                             MatrixMarketInfo* info = nullptr);

// Write in coordinate real format.  When symmetric_storage is true only the
// lower triangle is written and the header says "symmetric"; the caller is
// responsible for the matrix actually being symmetric.  Values are written
// with 17 significant digits so the round trip is exact.
void write_matrix_market(std::ostream& out, const CscMatrix& a,
                         bool symmetric_storage = false,
                         const std::string& comment = "");
void write_matrix_market(const std::string& path, const CscMatrix& a,
                         bool symmetric_storage = false,
                         const std::string& comment = "");

std::string to_string(MatrixMarketInfo::Field f);
std::string to_string(MatrixMarketInfo::Symmetry s);

