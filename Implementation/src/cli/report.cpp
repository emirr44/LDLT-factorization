#include "report.hpp"

#include <chrono>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "csc.hpp"
#include "matrix_market.hpp"
#include "numeric.hpp"
#include "refine.hpp"
#include "symbolic.hpp"
#include "validate.hpp"

namespace {

constexpr Real kStaticPivot = 1e-12;
constexpr Index kRefineSteps = 10;

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count() * 1e3;
}

std::string base_name(const std::string& path) {
    std::string b = std::filesystem::path(path).filename().string();
    const auto d = b.rfind(".mtx");
    if (d != std::string::npos && d + 4 == b.size()) b.erase(d);
    return b;
}

void row(std::ostringstream& os, const char* label, const std::string& value) {
    os << std::left << std::setw(28) << label << value << '\n';
}

std::string fixed2(double v) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(2) << v;
    return os.str();
}

std::string sci(double v, int prec = 3) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

}

FactorResult factorize_matrix(const std::string& path, Ordering ordering) {
    FactorResult res;
    res.matrix = base_name(path);
    res.suggested_file_name = res.matrix + "_" + to_string(ordering) + ".txt";

    std::ostringstream s;
    try {
        auto t0 = Clock::now();
        MatrixMarketInfo info;
        const CscMatrix a = read_matrix_market(path, true, &info);
        const double t_read = ms_since(t0);

        if (!a.is_square()) throw InvalidMatrix("matrix is not square");
        const bool symmetric = a.is_symmetric(1e-12);

        t0 = Clock::now();
        std::vector<Index> perm;
        if (ordering != Ordering::Natural) perm = compute_ordering(a, ordering);
        const double t_order = ms_since(t0);

        t0 = Clock::now();
        const Symbolic sym = Symbolic::analyze_fast(a, perm);
        const double t_sym = ms_since(t0);

        NumericOptions nopt;
        nopt.static_pivot_relative = kStaticPivot;

        t0 = Clock::now();
        const Numeric num = Numeric::factorize(a, sym, nopt);
        const double t_num = ms_since(t0);

        res.ok = num.ok();

        row(s, "matrix", res.matrix);
        row(s, "file", path);
        row(s, "header", "coordinate " + to_string(info.field) + ' ' +
                             to_string(info.symmetry));
        row(s, "n", std::to_string(a.nrows()));
        row(s, "nnz(A) stored", std::to_string(a.nnz()));
        row(s, "ordering", to_string(ordering));
        row(s, "nnz(L)", std::to_string(sym.nnz_L()));
        row(s, "fill nnz(L+D)/nnz(triu A)", fixed2(sym.fill_ratio(a)));
        row(s, "flops (numeric)", sci(sym.flops(), 4));
        s << '\n';
        row(s, "read", fixed2(t_read) + " ms");
        row(s, "ordering", fixed2(t_order) + " ms");
        row(s, "symbolic", fixed2(t_sym) + " ms");
        row(s, "numeric", fixed2(t_num) + " ms");

        if (!res.ok) {
            s << '\n';
            row(s, "status", "FAILED - zero pivot at k = " +
                                 std::to_string(num.failed_pivot()));
        }
        else {
            t0 = Clock::now();
            RefineOptions ro;
            ro.max_iterations = kRefineSteps;
            const RefinedSolveCheck rc = check_solve_refined(a, num, ro, true);
            const double t_solve = ms_since(t0);

            row(s, "solve + refine", fixed2(t_solve) + " ms");
            s << '\n';
            row(s, "status", "ok");
            row(s, "pivots", std::to_string(num.num_positive_pivots()) + " positive, " +
                                 std::to_string(num.num_negative_pivots()) + " negative");
            row(s, "min |d| / max |d|",
                sci(num.min_abs_pivot()) + " / " + sci(num.max_abs_pivot()));
            row(s, "static pivoting",
                "threshold " + sci(num.pivot_threshold()) + ", " +
                    std::to_string(num.num_perturbed_pivots()) + " perturbed" +
                    (num.num_perturbed_pivots() > 0
                         ? ", max " + sci(num.max_perturbation())
                         : std::string()));
            row(s, "backward error", sci(rc.before.residual.backward_error()));
            const Index steps = rc.refine.iterations;
            row(s, "  after refinement",
                sci(rc.after.residual.backward_error()) + "  (" +
                    std::to_string(steps) + (steps == 1 ? " step, " : " steps, ") +
                    (rc.refine.converged ? "converged" : "not converged") + ")");
            row(s, "forward error", sci(rc.after.forward_error));
        }

        if (!symmetric) {
            s << '\n' << "warning: matrix is not numerically symmetric;"
                         " only the upper triangle was used\n";
        }

        res.summary = s.str();

        std::ostringstream r;
        r << res.summary << '\n';
        if (res.ok) {
            r << "L (unit diagonal not stored):\n";
            num.L().print(r, 20, std::numeric_limits<Index>::max());
            r << "\nD:\n";
            r << std::setprecision(17);
            const std::vector<Real>& d = num.D();
            for (std::size_t i = 0; i < d.size(); ++i) {
                r << std::setw(10) << i << "  " << d[i] << '\n';
            }
            if (!perm.empty()) {
                r << "\nP (P[k] = j means column j of A becomes pivot k):\n";
                for (std::size_t k = 0; k < perm.size(); ++k) {
                    r << std::setw(10) << k << "  " << perm[k] << '\n';
                }
            }
        }
        else {
            r << "no factors: the factorization stopped at a zero pivot\n";
        }
        res.report = r.str();
    }
    catch (const std::exception& e) {
        res.ok = false;
        res.threw = true;
        res.error = e.what();
        std::ostringstream err;
        row(err, "matrix", res.matrix);
        row(err, "file", path);
        err << '\n';
        row(err, "status", std::string("ERROR - ") + e.what());
        res.summary = err.str();
        res.report = res.summary;
    }

    return res;
}
