// ldlt/etree.hpp - elimination tree, postorder, column counts, symmetric
// permutation.  Ports of cs_etree, cs_post, cs_counts, cs_symperm from
// CSparse (Davis, "Direct Methods for Sparse Linear Systems", ch. 4), LGPL.
//
// These replace the O(nnz(L)) path-following symbolic analysis of ldl.c with
// algorithms that run in time nearly linear in nnz(A):
//   elimination_tree   O(nnz(A) * alpha)  via path compression (ancestor array)
//   postorder          O(n)
//   column_counts      O(nnz(A) * alpha)  Gilbert/Ng/Peyton skeleton algorithm
//
// All routines read only the upper triangular part (i <= j) of the matrix, so
// the full symmetric matrix or its upper triangle may be passed.  Columns need
// not be sorted.
#pragma once

#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/types.hpp"


// parent[j] = smallest i > j with L(i,j) != 0, or -1 if j is a root.
std::vector<Index> elimination_tree(const CscMatrix& a);

// Postorder of the forest given by parent[]: post[k] is the kth node visited.
// Children are visited before their parent; siblings in increasing order.
std::vector<Index> postorder(const std::vector<Index>& parent);

// Number of nonzeros in each column of L, EXCLUDING the diagonal (matches
// Symbolic::col_counts and ldl.c's Lnz).  Requires the etree and a postorder.
std::vector<Index> column_counts(const CscMatrix& a, const std::vector<Index>& parent,
                                 const std::vector<Index>& post);

// C = P A P' restricted to the upper triangle, where pinv is the inverse
// permutation (pinv[j] = k means row/col j of A becomes row/col k of C).
// Only entries with i <= j of A are read.  Columns of C are not sorted.
CscMatrix permute_symmetric_upper(const CscMatrix& a, const std::vector<Index>& pinv);

// pinv[p[k]] = k
std::vector<Index> invert_perm(const std::vector<Index>& p);

// q[k] = p[post[k]]: the ordering "p followed by post".  If p is empty
// (identity) the result is post itself.
std::vector<Index> compose_perm(const std::vector<Index>& p, const std::vector<Index>& post);

