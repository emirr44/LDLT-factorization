// types.hpp - basic scalar types and error classes shared by the library.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>


// index type (for all dimensions, pointers and row indices), 64-bit so that nnz(L) cannot overflow on large factorizations
using Index = std::int64_t;

// scalar type of matrix entries.
using Real = double;

class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatrix : public Error {
public:
    using Error::Error;
};

class IoError : public Error {
public:
    using Error::Error;
};

