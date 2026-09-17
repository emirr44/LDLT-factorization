// report.hpp - runs the LDL' pipeline on one file and formats the results.
//
// Pure computation: nothing is printed here and no exception escapes, so the
// caller decides what to show and where to write it.
#pragma once

#include <string>

#include "ldlt/ordering.hpp"


struct FactorResult {
    bool ok = false;           // factorization completed (no zero pivot)
    bool threw = false;        // read / analysis error, see `error`
    std::string matrix;        // matrix name
    std::string error;         // exception message when threw
    std::string summary;       // the short report shown on screen
    std::string report;        // summary + the complete L, D and permutation
    std::string suggested_file_name;   // e.g. "bcsstk14_amd.txt"
};

// Static pivoting and refinement are always on: 1e-12 was safe across the whole
// test set and is a no-op on well-scaled SPD matrices. (1e-8 is the classic
// sqrt(eps) choice but degrades badly scaled ones such as bcsstk13/14.)
FactorResult factorize_matrix(const std::string& path, Ordering ordering);
