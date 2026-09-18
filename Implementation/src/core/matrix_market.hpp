// Matrix Market (.mtx) reader and writer.
#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "coo.hpp"
#include "csc.hpp"


struct MatrixMarketInfo {
    enum class Field { Double, Integer, Pattern };
    enum class Symmetry { General, Symmetric, SkewSymmetric };

    Field field = Field::Double;
    Symmetry symmetry = Symmetry::General;
    Index nrows = 0;
    Index ncols = 0;
    Index nnz_in_file = 0;
    std::vector<std::string> comments;

    bool is_symmetric_storage() const { return symmetry != Symmetry::General; }
};

CooMatrix read_matrix_market_coo(std::istream& in,
                                 bool expand_symmetric = true,
                                 MatrixMarketInfo* info = nullptr);
CooMatrix read_matrix_market_coo(const std::string& path,
                                 bool expand_symmetric = true,
                                 MatrixMarketInfo* info = nullptr);

CscMatrix read_matrix_market(std::istream& in,
                             bool expand_symmetric = true,
                             MatrixMarketInfo* info = nullptr);
CscMatrix read_matrix_market(const std::string& path,
                             bool expand_symmetric = true,
                             MatrixMarketInfo* info = nullptr);

void write_matrix_market(std::ostream& out, const CscMatrix& a,
                         bool symmetric_storage = false,
                         const std::string& comment = "");
void write_matrix_market(const std::string& path, const CscMatrix& a,
                         bool symmetric_storage = false,
                         const std::string& comment = "");

std::string to_string(MatrixMarketInfo::Field f);
std::string to_string(MatrixMarketInfo::Symmetry s);

