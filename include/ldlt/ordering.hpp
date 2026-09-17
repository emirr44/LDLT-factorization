// ldlt/ordering.hpp - fill-reducing orderings for symmetric matrices.
//
// All orderings return a permutation P in the convention used everywhere in
// this library: P[k] = j means row/column j of A becomes pivot k, so the
// factorization is P A P' = L D L'.  Only the pattern of A is used; the full
// symmetric matrix or either triangle may be passed (the pattern is
// symmetrized internally and the diagonal is ignored).
//
//   Natural        identity
//   RCM            reverse Cuthill-McKee: bandwidth/profile reduction.  BFS
//                  from a pseudo-peripheral node per component, neighbours
//                  visited by increasing degree, order reversed.  Cheap
//                  baseline; good for banded-ish problems, poor in general.
//   MinimumDegree  exact minimum degree on an explicit elimination graph.
//                  Simple and easy to trust, but fill edges are materialized,
//                  so time and memory grow with nnz(L).  Reference for AMD.
//   AMD            approximate minimum degree on a quotient graph (Amestoy,
//                  Davis, Duff 1996): elements instead of fill edges, element
//                  absorption, supervariable detection, and the approximate
//                  external degree bound in place of the exact degree.  This
//                  is a readable re-implementation, not SuiteSparse AMD; it
//                  omits dense-row detection and aggressive absorption.
#pragma once

#include <string>
#include <vector>

#include "ldlt/csc.hpp"
#include "ldlt/types.hpp"


enum class Ordering { Natural, RCM, MinimumDegree, AMD };

const char* to_string(Ordering o);
// Accepts "natural", "rcm", "md", "amd" (case-insensitive).  False if unknown.
bool parse_ordering(const std::string& s, Ordering& out);

std::vector<Index> compute_ordering(const CscMatrix& a, Ordering method);

std::vector<Index> natural_ordering(Index n);
std::vector<Index> rcm_ordering(const CscMatrix& a);
std::vector<Index> minimum_degree_ordering(const CscMatrix& a);
std::vector<Index> amd_ordering(const CscMatrix& a);

// Adjacency structure of the symmetrized pattern of A without the diagonal.
// Neighbours of node i are adj[ptr[i] .. ptr[i+1]-1], sorted, no duplicates.
struct Graph {
    Index n = 0;
    std::vector<Index> ptr;
    std::vector<Index> adj;

    static Graph from_matrix(const CscMatrix& a);
    Index degree(Index i) const { return ptr[static_cast<std::size_t>(i) + 1] - ptr[static_cast<std::size_t>(i)]; }
    Index num_edges() const { return static_cast<Index>(adj.size()) / 2; }
};

// Fill statistics of an ordering, from the symbolic analysis alone.
struct FillStats {
    Index nnz_L = 0;        // strictly below the diagonal
    double flops = 0.0;     // numeric-phase flop count
    double fill_ratio = 0.0;
};
FillStats fill_for_ordering(const CscMatrix& a, const std::vector<Index>& perm);

