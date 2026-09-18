// fill-reducing orderings for symmetric matrices.

/*
Natural- identity
RCM - reverse Cuthill-McKee: bandwidth/profile reduction
MinimumDegree  - exact minimum degree on an explicit elimination graph
AMD - approximate minimum degree on a quotient graph (Amestoy, Davis, Duff 1996): elements instead of fill edges, element
absorption, supervariable detection, and the approximate external degree bound in place of the exact degree
*/


#pragma once

#include <string>
#include <vector>

#include "csc.hpp"
#include "types.hpp"


enum class Ordering { Natural, RCM, MinimumDegree, AMD };

const char* to_string(Ordering o);
bool parse_ordering(const std::string& s, Ordering& out);

std::vector<Index> compute_ordering(const CscMatrix& a, Ordering method);

std::vector<Index> natural_ordering(Index n);
std::vector<Index> rcm_ordering(const CscMatrix& a);
std::vector<Index> minimum_degree_ordering(const CscMatrix& a);
std::vector<Index> amd_ordering(const CscMatrix& a);

struct Graph {
    Index n = 0;
    std::vector<Index> ptr;
    std::vector<Index> adj;

    static Graph from_matrix(const CscMatrix& a);
    Index degree(Index i) const { return ptr[static_cast<std::size_t>(i) + 1] - ptr[static_cast<std::size_t>(i)]; }
    Index num_edges() const { return static_cast<Index>(adj.size()) / 2; }
};

struct FillStats {
    Index nnz_L = 0;
    double flops = 0.0;     
    double fill_ratio = 0.0;
};
FillStats fill_for_ordering(const CscMatrix& a, const std::vector<Index>& perm);

