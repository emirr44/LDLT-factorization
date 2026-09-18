// elimination tree, postorder, column counts, symmetric permutation

// ports of cs_etree, cs_post, cs_counts, cs_symperm from
// CSparse (Davis, "Direct Methods for Sparse Linear Systems", ch. 4), LGPL.

#pragma once

#include <vector>

#include "csc.hpp"
#include "types.hpp"


std::vector<Index> elimination_tree(const CscMatrix& a);

// postorder of the forest given by parent[]: post[k] is the kth node visited
std::vector<Index> postorder(const std::vector<Index>& parent);

// number of nonzeros in each column of L, EXCLUDING the diagonal
std::vector<Index> column_counts(const CscMatrix& a, const std::vector<Index>& parent,
                                 const std::vector<Index>& post);

// C = P A P' restricted to the upper triangle
CscMatrix permute_symmetric_upper(const CscMatrix& a, const std::vector<Index>& pinv);

// pinv[p[k]] = k
std::vector<Index> invert_perm(const std::vector<Index>& p);

// q[k] = p[post[k]]: the ordering "p followed by post"
std::vector<Index> compose_perm(const std::vector<Index>& p, const std::vector<Index>& post);

