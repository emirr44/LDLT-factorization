// ldlt/types.hpp - basic scalar types and error classes shared by the library.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>


// Index type used for all dimensions, pointers and row indices.
// 64-bit so that nnz(L) cannot overflow on large factorizations.
using Index = std::int64_t;

// Scalar type of matrix entries.
using Real = double;

// Base class for every exception thrown by this library.
class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Thrown when a matrix (COO or CSC) violates its structural invariants.
class InvalidMatrix : public Error {
public:
    using Error::Error;
};

// Thrown on file / parse problems (e.g. Matrix Market reader).
class IoError : public Error {
public:
    using Error::Error;
};

