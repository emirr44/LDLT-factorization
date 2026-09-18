#include "ordering.hpp"

#include <algorithm>
#include <cctype>
#include <numeric>
#include <queue>
#include <utility>

#include "coo.hpp"
#include "symbolic.hpp"


namespace {
std::size_t sz(Index n) { return static_cast<std::size_t>(n); }
}

// -- names --

const char* to_string(Ordering o) {
    switch (o) {
        case Ordering::Natural: return "natural";
        case Ordering::RCM: return "rcm";
        case Ordering::MinimumDegree: return "md";
        case Ordering::AMD: return "amd";
    }
    return "?";
}

bool parse_ordering(const std::string& s, Ordering& out) {
    std::string t;
    for (char c : s) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "natural" || t == "none" || t == "identity") { out = Ordering::Natural; return true; }
    if (t == "rcm") { out = Ordering::RCM; return true; }
    if (t == "md" || t == "mindeg" || t == "minimum-degree") { out = Ordering::MinimumDegree; return true; }
    if (t == "amd") { out = Ordering::AMD; return true; }
    return false;
}

std::vector<Index> compute_ordering(const CscMatrix& a, Ordering method) {
    switch (method) {
        case Ordering::Natural: return natural_ordering(a.nrows());
        case Ordering::RCM: return rcm_ordering(a);
        case Ordering::MinimumDegree: return minimum_degree_ordering(a);
        case Ordering::AMD: return amd_ordering(a);
    }
    throw Error("compute_ordering: unknown method");
}

std::vector<Index> natural_ordering(Index n) {
    std::vector<Index> p(sz(n));
    std::iota(p.begin(), p.end(), 0);
    return p;
}

FillStats fill_for_ordering(const CscMatrix& a, const std::vector<Index>& perm) {
    const Symbolic s = Symbolic::analyze_fast(a, perm);
    FillStats f;
    f.nnz_L = s.nnz_L();
    f.flops = s.flops();
    f.fill_ratio = s.fill_ratio(a);
    return f;
}

// -- graph --

Graph Graph::from_matrix(const CscMatrix& a) {
    if (!a.is_square()) throw InvalidMatrix("Graph::from_matrix: matrix must be square");
    const Index n = a.nrows();
    // symmetrize the pattern via COO -> CSC, which sorts and removes duplicates
    CooMatrix coo(n, n);
    coo.reserve(2 * a.nnz());
    const auto& Ap = a.colptr();
    const auto& Ai = a.rowind();
    for (Index j = 0; j < n; ++j) {
        for (Index p = Ap[sz(j)]; p < Ap[sz(j) + 1]; ++p) {
            const Index i = Ai[sz(p)];
            if (i == j) continue;
            coo.add(i, j, 1.0);
            coo.add(j, i, 1.0);
        }
    }
    const CscMatrix s = CscMatrix::from_coo(coo);
    Graph g;
    g.n = n;
    g.ptr = s.colptr();
    g.adj = s.rowind();
    return g;
}

// -- reverse Cuthill-McKee --
namespace {

/*
breadth-first search from root over nodes with level[] == -1.  Appends the
visited nodes to `visited` in BFS order, sets their level, and returns the
eccentricity (largest level)
*/ 

Index bfs_levels(const Graph& g, Index root, std::vector<Index>& level, std::vector<Index>& visited) {
    const std::size_t start = visited.size();
    level[sz(root)] = 0;
    visited.push_back(root);
    Index ecc = 0;
    for (std::size_t q = start; q < visited.size(); ++q) {
        const Index v = visited[q];
        for (Index p = g.ptr[sz(v)]; p < g.ptr[sz(v) + 1]; ++p) {
            const Index u = g.adj[sz(p)];
            if (level[sz(u)] != -1) continue;
            level[sz(u)] = level[sz(v)] + 1;
            ecc = std::max(ecc, level[sz(u)]);
            visited.push_back(u);
        }
    }
    return ecc;
}

// George-Liu pseudo-peripheral node of the component containing `start`.
Index pseudo_peripheral(const Graph& g, Index start, std::vector<Index>& level) {
    std::vector<Index> visited;
    Index x = start;
    Index ecc_x = bfs_levels(g, x, level, visited);
    for (int iter = 0; iter < 10; ++iter) {
        // candidate: minimum-degree node in the last level
        Index y = -1;
        for (Index v : visited) {
            if (level[sz(v)] != ecc_x) continue;
            if (y == -1 || g.degree(v) < g.degree(y) || (g.degree(v) == g.degree(y) && v < y)) y = v;
        }
        for (Index v : visited) level[sz(v)] = -1;   // reset for the next BFS
        visited.clear();
        const Index ecc_y = bfs_levels(g, y, level, visited);
        if (ecc_y > ecc_x) {
            x = y;
            ecc_x = ecc_y;
        } else {
            break;
        }
    }
    for (Index v : visited) level[sz(v)] = -1;
    return x;
}

}  // namespace

std::vector<Index> rcm_ordering(const CscMatrix& a) {
    const Graph g = Graph::from_matrix(a);
    const Index n = g.n;
    std::vector<Index> order;
    order.reserve(sz(n));
    std::vector<Index> level(sz(n), -1);    // scratch for pseudo-peripheral search
    std::vector<char> visited(sz(n), 0);
    std::vector<Index> nbrs;

    for (Index s = 0; s < n; ++s) {
        if (visited[sz(s)]) continue;
        const Index root = pseudo_peripheral(g, s, level);
        // Cuthill-McKee BFS from root, neighbours by increasing degree
        const std::size_t head = order.size();
        visited[sz(root)] = 1;
        order.push_back(root);
        for (std::size_t q = head; q < order.size(); ++q) {
            const Index v = order[q];
            nbrs.clear();
            for (Index p = g.ptr[sz(v)]; p < g.ptr[sz(v) + 1]; ++p) {
                const Index u = g.adj[sz(p)];
                if (!visited[sz(u)]) {
                    visited[sz(u)] = 1;
                    nbrs.push_back(u);
                }
            }
            std::sort(nbrs.begin(), nbrs.end(), [&](Index x, Index y) {
                const Index dx = g.degree(x), dy = g.degree(y);
                return dx != dy ? dx < dy : x < y;
            });
            order.insert(order.end(), nbrs.begin(), nbrs.end());
        }
    }
    std::reverse(order.begin(), order.end());
    return order;
}

// -- exact minimum degree (elimination graph) --

std::vector<Index> minimum_degree_ordering(const CscMatrix& a) {
    const Graph g = Graph::from_matrix(a);
    const Index n = g.n;

    std::vector<std::vector<Index>> adj(sz(n));
    std::vector<Index> degree(sz(n));
    for (Index i = 0; i < n; ++i) {
        adj[sz(i)].assign(g.adj.begin() + g.ptr[sz(i)], g.adj.begin() + g.ptr[sz(i) + 1]);
        degree[sz(i)] = static_cast<Index>(adj[sz(i)].size());
    }
    std::vector<char> eliminated(sz(n), 0);

    // lazy min-heap of (degree, node); stale entries are skipped on pop
    using Entry = std::pair<Index, Index>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> heap;
    for (Index i = 0; i < n; ++i) heap.emplace(degree[sz(i)], i);

    std::vector<Index> order;
    order.reserve(sz(n));
    std::vector<Index> merged;

    while (static_cast<Index>(order.size()) < n) {
        const Entry top = heap.top();
        heap.pop();
        const Index v = top.second;
        if (eliminated[sz(v)] || top.first != degree[sz(v)]) continue;   // stale
        eliminated[sz(v)] = 1;
        order.push_back(v);

        std::vector<Index> nb = std::move(adj[sz(v)]);   // all non-eliminated by invariant
        adj[sz(v)] = std::vector<Index>();
        for (Index u : nb) {
            // adj[u] = (adj[u] \ {v}) U (nb \ {u}), both sorted
            merged.clear();
            const std::vector<Index>& au = adj[sz(u)];
            std::size_t x = 0, y = 0;
            while (x < au.size() || y < nb.size()) {
                if (x < au.size() && au[x] == v) { ++x; continue; }
                if (y < nb.size() && nb[y] == u) { ++y; continue; }
                if (y == nb.size() || (x < au.size() && au[x] < nb[y])) merged.push_back(au[x++]);
                else if (x == au.size() || nb[y] < au[x]) merged.push_back(nb[y++]);
                else { merged.push_back(au[x]); ++x; ++y; }   // equal: keep one
            }
            adj[sz(u)].swap(merged);
            degree[sz(u)] = static_cast<Index>(adj[sz(u)].size());
            heap.emplace(degree[sz(u)], u);
        }
    }
    return order;
}

// approximate minimum degree (quotient graph)

/*
Quotient graph: each not-yet-eliminated node is a *variable* i with
   A_i : adjacent variables through original (still uncovered) edges
   E_i : adjacent elements
 and each eliminated pivot p is an *element* with
  L_p : the variables adjacent to p when it was eliminated (= the nonzero
        pattern of column p of L, as supervariables)
*/

std::vector<Index> amd_ordering(const CscMatrix& a) {
    const Graph g = Graph::from_matrix(a);
    const Index n = g.n;
    if (n == 0) return {};

    std::vector<std::vector<Index>> A(sz(n)), E(sz(n)), L(sz(n));
    std::vector<std::vector<Index>> members(sz(n));   // variables merged into i
    std::vector<Index> nv(sz(n), 1);                  // supervariable weight (0 = merged away)
    std::vector<Index> degree(sz(n));                 // approximate external degree
    std::vector<Index> esize(sz(n), 0);               // weighted |L_e|
    enum : char { VARIABLE = 0, ELEMENT = 1, MERGED = 2 };
    std::vector<char> status(sz(n), VARIABLE);
    std::vector<char> alive(sz(n), 0);                // element not yet absorbed

    // degree buckets (doubly linked lists), all principal variables are in one
    std::vector<Index> head(sz(n), -1), nxt(sz(n), -1), prv(sz(n), -1);
    auto bucket_insert = [&](Index i, Index d) {
        nxt[sz(i)] = head[sz(d)];
        prv[sz(i)] = -1;
        if (head[sz(d)] != -1) prv[sz(head[sz(d)])] = i;
        head[sz(d)] = i;
    };
    auto bucket_remove = [&](Index i, Index d) {
        if (prv[sz(i)] != -1) nxt[sz(prv[sz(i)])] = nxt[sz(i)];
        else head[sz(d)] = nxt[sz(i)];
        if (nxt[sz(i)] != -1) prv[sz(nxt[sz(i)])] = prv[sz(i)];
    };

    for (Index i = 0; i < n; ++i) {
        A[sz(i)].assign(g.adj.begin() + g.ptr[sz(i)], g.adj.begin() + g.ptr[sz(i) + 1]);
        degree[sz(i)] = static_cast<Index>(A[sz(i)].size());
        bucket_insert(i, degree[sz(i)]);
    }

    std::vector<Index> mark(sz(n), 0);
    Index stamp = 0;
    std::vector<Index> w(sz(n), 0), wstamp(sz(n), 0);   // w(e) = |L_e \ L_p|
    Index wflg = 0;
    std::vector<Index> hash(sz(n), 0), hhead(sz(n), -1), hnext(sz(n), -1);
    std::vector<Index> Lp;
    std::vector<Index> order;
    order.reserve(sz(n));

    Index k = 0;
    Index mindeg = 0;
    while (k < n) {
        // - select pivot p of minimum approximate degree -
        while (mindeg < n && head[sz(mindeg)] == -1) ++mindeg;
        if (mindeg >= n) throw Error("amd_ordering: internal error, no variable left");
        const Index p = head[sz(mindeg)];
        bucket_remove(p, mindeg);

        ++stamp;
        mark[sz(p)] = stamp;
        Lp.clear();
        Index lp_weight = 0;
        auto add_to_Lp = [&](Index i) {
            if (status[sz(i)] != VARIABLE || nv[sz(i)] == 0 || mark[sz(i)] == stamp) return;
            mark[sz(i)] = stamp;
            Lp.push_back(i);
            lp_weight += nv[sz(i)];
        };
        for (Index i : A[sz(p)]) add_to_Lp(i);
        for (Index e : E[sz(p)]) {
            if (!alive[sz(e)]) continue;
            for (Index i : L[sz(e)]) add_to_Lp(i);
            alive[sz(e)] = 0;               
            std::vector<Index>().swap(L[sz(e)]);
        }
        std::vector<Index>().swap(A[sz(p)]);
        std::vector<Index>().swap(E[sz(p)]);

        status[sz(p)] = ELEMENT;
        alive[sz(p)] = 1;
        esize[sz(p)] = lp_weight;
        order.push_back(p);
        for (Index j : members[sz(p)]) order.push_back(j);
        std::vector<Index>().swap(members[sz(p)]);
        k += nv[sz(p)];
        nv[sz(p)] = 0;

        // - w(e) = |L_e \ L_p| for every element adjacent to L_p -
        ++wflg;
        for (Index i : Lp) {
            for (Index e : E[sz(i)]) {
                if (!alive[sz(e)]) continue;
                if (wstamp[sz(e)] != wflg) {
                    wstamp[sz(e)] = wflg;
                    w[sz(e)] = esize[sz(e)];
                }
                w[sz(e)] -= nv[sz(i)];
            }
        }

        // - update each i in L_p: prune lists, approximate degree, hash -
        for (Index i : Lp) {
            bucket_remove(i, degree[sz(i)]);
            std::vector<Index>& Ei = E[sz(i)];
            std::vector<Index>& Ai = A[sz(i)];

            Index deg_e = 0;
            Index h = 0;
            std::size_t out = 0;
            for (std::size_t q = 0; q < Ei.size(); ++q) {
                const Index e = Ei[q];
                if (!alive[sz(e)]) continue;
                Ei[out++] = e;
                deg_e += w[sz(e)];
                h += e;
            }
            Ei.resize(out);
            Ei.push_back(p);
            h += p;

            Index a_weight = 0;
            out = 0;
            for (std::size_t q = 0; q < Ai.size(); ++q) {
                const Index j = Ai[q];
                if (status[sz(j)] != VARIABLE || nv[sz(j)] == 0 || mark[sz(j)] == stamp) continue;
                Ai[out++] = j;
                a_weight += nv[sz(j)];
                h += j;
            }
            Ai.resize(out);

            const Index ext_lp = lp_weight - nv[sz(i)];
            Index d = a_weight + ext_lp + deg_e;
            d = std::min(d, degree[sz(i)] + ext_lp);
            d = std::min(d, n - k - nv[sz(i)]);
            degree[sz(i)] = std::max<Index>(d, 0);

            std::sort(Ei.begin(), Ei.end());
            std::sort(Ai.begin(), Ai.end());
            hash[sz(i)] = ((h % n) + n) % n;
        }

        for (Index i : Lp) {
            const Index h = hash[sz(i)];
            hnext[sz(i)] = hhead[sz(h)];
            hhead[sz(h)] = i;
        }
        for (Index i : Lp) {
            const Index h = hash[sz(i)];
            if (hhead[sz(h)] == -1) continue;
            for (Index x = hhead[sz(h)]; x != -1; x = hnext[sz(x)]) {
                if (nv[sz(x)] == 0) continue;
                for (Index y = hnext[sz(x)]; y != -1; y = hnext[sz(y)]) {
                    if (nv[sz(y)] == 0) continue;
                    if (E[sz(x)] != E[sz(y)] || A[sz(x)] != A[sz(y)]) continue;
                    degree[sz(x)] -= nv[sz(y)];
                    nv[sz(x)] += nv[sz(y)];
                    nv[sz(y)] = 0;
                    status[sz(y)] = MERGED;
                    members[sz(x)].push_back(y);
                    members[sz(x)].insert(members[sz(x)].end(), members[sz(y)].begin(), members[sz(y)].end());
                    std::vector<Index>().swap(members[sz(y)]);
                    std::vector<Index>().swap(A[sz(y)]);
                    std::vector<Index>().swap(E[sz(y)]);
                }
            }
            hhead[sz(h)] = -1;
        }

        // - finalize: L_p keeps principal variables, re-insert degrees -
        std::vector<Index>& lp_final = L[sz(p)];
        lp_final.clear();
        for (Index i : Lp) {
            if (nv[sz(i)] == 0) continue;
            lp_final.push_back(i);
            bucket_insert(i, degree[sz(i)]);
            mindeg = std::min(mindeg, degree[sz(i)]);
        }
        if (lp_final.empty()) alive[sz(p)] = 0;
    }

    if (static_cast<Index>(order.size()) != n) throw Error("amd_ordering: internal error, incomplete ordering");
    return order;
}

