
#include "batch.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "csc.hpp"
#include "matrix_market.hpp"
#include "numeric.hpp"
#include "ordering.hpp"
#include "refine.hpp"
#include "symbolic.hpp"
#include "validate.hpp"

namespace {

void print_usage(std::ostream& os) {
    os <<
        "usage:\n"
        "  ldlt interactive session (start here)\n"
        "\n"
        "  ldlt info    <file.mtx>\n"
        "  ldlt print   <file.mtx> [dense_limit]\n"
        "  ldlt convert <in.mtx> <out.mtx> [--symmetric]\n"
        "  ldlt factor  <file.mtx>... [--order natural|rcm|md|amd] [--postorder]\n"
        "               [--print] [--check-factor] [--ones] [--tol <t>]\n"
        "               [--static-pivot <rel>] [--refine <max_iters>]\n"
        "  ldlt orderings <file.mtx>...\n"
        "\n"
        "options for every mode:\n"
        "  --plain        no colour or box drawing\n"
        "  --root <dir>   project folder holding data/ (default: where ldlt was built)\n";
}

int usage() {
    print_usage(std::cerr);
    return 2;
}

void print_header(const std::string& path, const MatrixMarketInfo& info) {
    std::cout << "file                  : " << path << '\n'
              << "header                : coordinate " << to_string(info.field) << ' '
              << to_string(info.symmetry) << '\n'
              << "entries in file       : " << info.nnz_in_file << '\n';
}

int cmd_info(const std::vector<std::string>& args) {
    if (args.size() != 1) return usage();
    MatrixMarketInfo info;
    const CscMatrix a = read_matrix_market(args[0], true, &info);
    print_header(args[0], info);
    std::cout << a.stats();
    if (a.is_square() && !a.is_symmetric(1e-12)) {
        std::cout << "note: matrix is not symmetric; LDL' factorization needs a symmetric matrix\n";
    }
    return 0;
}

int cmd_print(const std::vector<std::string>& args) {
    if (args.empty() || args.size() > 2) return usage();
    Index limit = 20;
    if (args.size() == 2) limit = std::stoll(args[1]);
    MatrixMarketInfo info;
    const CscMatrix a = read_matrix_market(args[0], true, &info);
    print_header(args[0], info);
    a.print(std::cout, limit);
    return 0;
}

int cmd_convert(const std::vector<std::string>& args) {
    if (args.size() < 2 || args.size() > 3) return usage();
    bool symmetric = false;
    if (args.size() == 3) {
        if (args[2] != "--symmetric") return usage();
        symmetric = true;
    }
    const CscMatrix a = read_matrix_market(args[0]);
    if (symmetric && !a.is_symmetric(1e-12)) {
        std::cerr << "error: matrix is not symmetric, refusing to write symmetric storage\n";
        return 1;
    }
    write_matrix_market(args[1], a, symmetric, "converted by ldlt_cli from " + args[0]);
    std::cout << "wrote " << args[1] << " (" << a.nrows() << " x " << a.ncols() << ", "
              << (symmetric ? a.lower_triangle(true).nnz() : a.nnz()) << " entries)\n";
    return 0;
}

// -- factor --

struct FactorOptions {
    Ordering order = Ordering::Natural;
    bool postorder = false;
    bool print = false;
    bool check_factor = false;
    bool ones = false;
    Real tol = 0.0;
    Real static_pivot = 0.0;   // relative threshold, 0 = off
    Index refine = 0;          // max refinement steps, 0 = off
};

struct FactorRun {
    std::string name;
    Index n = 0, nnz_a = 0, nnz_l = 0;
    double fill = 0.0, flops = 0.0;
    double t_read = 0.0, t_order = 0.0, t_sym = 0.0, t_num = 0.0, t_solve = 0.0;
    bool ok = false;
    Index failed_pivot = -1, neg = 0, perturbed = 0;
    Real dmin = 0.0, dmax = 0.0;
    Real backward = 0.0, forward = 0.0, recon = -1.0;
    Real backward0 = 0.0;
    Index refine_iters = 0;
    bool refine_converged = false;
};

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

std::string base_name(const std::string& path) {
    const auto p = path.find_last_of("/\\");
    std::string b = (p == std::string::npos) ? path : path.substr(p + 1);
    const auto d = b.rfind(".mtx");
    if (d != std::string::npos && d + 4 == b.size()) b.erase(d);
    return b;
}

FactorRun factor_one(const std::string& path, const FactorOptions& opt, bool verbose) {
    FactorRun r;
    r.name = base_name(path);

    auto t0 = Clock::now();
    MatrixMarketInfo info;
    const CscMatrix a = read_matrix_market(path, true, &info);
    r.t_read = seconds_since(t0);
    if (!a.is_square()) throw InvalidMatrix("factor: matrix is not square");
    r.n = a.nrows();
    r.nnz_a = a.nnz();
    if (verbose) {
        print_header(path, info);
        if (!a.is_symmetric(1e-12)) {
            std::cout << "warning: matrix is not numerically symmetric; only the upper triangle is used\n";
        }
    }

    t0 = Clock::now();
    std::vector<Index> perm;
    if (opt.order != Ordering::Natural) perm = compute_ordering(a, opt.order);
    r.t_order = seconds_since(t0);

    t0 = Clock::now();
    Symbolic sym = Symbolic::analyze_fast(a, perm);
    if (opt.postorder) sym = Symbolic::analyze_fast(a, sym.postordered_perm());
    r.t_sym = seconds_since(t0);
    r.nnz_l = sym.nnz_L();
    r.fill = sym.fill_ratio(a);
    r.flops = sym.flops();

    NumericOptions nopt;
    nopt.zero_pivot_tolerance = opt.tol;
    nopt.static_pivot_relative = opt.static_pivot;
    t0 = Clock::now();
    const Numeric num = Numeric::factorize(a, sym, nopt);
    r.t_num = seconds_since(t0);
    r.ok = num.ok();
    r.failed_pivot = num.failed_pivot();
    r.neg = num.num_negative_pivots();
    r.perturbed = num.num_perturbed_pivots();
    r.dmin = num.min_abs_pivot();
    r.dmax = num.max_abs_pivot();

    if (r.ok) {
        t0 = Clock::now();
        if (opt.refine > 0) {
            RefineOptions ro;
            ro.max_iterations = opt.refine;
            const RefinedSolveCheck rc = check_solve_refined(a, num, ro, !opt.ones);
            r.backward0 = rc.before.residual.backward_error();
            r.backward = rc.after.residual.backward_error();
            r.forward = rc.after.forward_error;
            r.refine_iters = rc.refine.iterations;
            r.refine_converged = rc.refine.converged;
        } else {
            const SolveCheck sc = check_solve(a, num, !opt.ones);
            r.backward0 = r.backward = sc.residual.backward_error();
            r.forward = sc.forward_error;
        }
        r.t_solve = seconds_since(t0);
        if (opt.check_factor && r.n <= 3000) r.recon = reconstruction_error(a, num);
    }

    if (verbose) {
        std::cout << std::left;
        std::cout << "n                     : " << r.n << '\n'
                  << "nnz(A) stored         : " << r.nnz_a << '\n'
                  << "ordering              : " << to_string(opt.order)
                  << (opt.postorder ? " + postorder" : "") << '\n'
                  << "nnz(L) below diagonal : " << r.nnz_l << '\n'
                  << "fill ratio nnz(L+D)/nnz(triu A): " << std::setprecision(3) << r.fill << '\n'
                  << "flops (numeric)       : " << std::setprecision(4) << r.flops << '\n'
                  << "time read / order / symbolic / numeric / solve+check [ms]: "
                  << std::fixed << std::setprecision(2)
                  << r.t_read * 1e3 << " / " << r.t_order * 1e3 << " / " << r.t_sym * 1e3
                  << " / " << r.t_num * 1e3 << " / " << r.t_solve * 1e3 << std::defaultfloat << '\n';
        if (!r.ok) {
            std::cout << "STATUS                : FAILED, zero pivot at k = " << r.failed_pivot << '\n';
        } else {
            std::cout << "status                : ok\n"
                      << "pivots                : " << (r.n - r.neg) << " positive, " << r.neg
                      << " negative  (min |d| = " << std::setprecision(4) << r.dmin
                      << ", max |d| = " << r.dmax << ")\n";
            if (opt.static_pivot > 0.0) {
                std::cout << "static pivoting       : threshold " << std::scientific << std::setprecision(3)
                          << num.pivot_threshold() << ", " << r.perturbed << " pivots perturbed"
                          << ", max perturbation " << num.max_perturbation() << std::defaultfloat << '\n';
            }
            std::cout << "backward error ||b-Ax|| / (||A|| ||x|| + ||b||) : "
                      << std::scientific << std::setprecision(3) << r.backward0 << '\n';
            if (opt.refine > 0) {
                std::cout << "  after refinement    : " << r.backward << "  (" << r.refine_iters
                          << " steps, " << (r.refine_converged ? "converged" : "not converged") << ")\n";
            }
            std::cout << "forward error  ||x-x_true|| / ||x_true||        : " << r.forward << '\n';
            if (r.recon >= 0.0) {
                std::cout << "factor error   ||PAP'-LDL'||_F / ||A||_F         : " << r.recon << '\n';
            }
            std::cout << std::defaultfloat;
        }
        if (opt.print || r.n <= 12) {
            std::cout << "\nL (unit diagonal not stored):\n";
            num.L().print(std::cout, opt.print ? 40 : 12);
            std::cout << "D:";
            std::cout << std::setprecision(6);
            for (Real d : num.D()) std::cout << ' ' << d;
            std::cout << '\n';
        }
    }
    return r;
}

int cmd_factor(const std::vector<std::string>& args) {
    FactorOptions opt;
    std::vector<std::string> files;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& s = args[i];
        if (s == "--print") opt.print = true;
        else if (s == "--check-factor") opt.check_factor = true;
        else if (s == "--ones") opt.ones = true;
        else if (s == "--postorder") opt.postorder = true;
        else if (s == "--order") {
            if (i + 1 >= args.size()) return usage();
            if (!parse_ordering(args[++i], opt.order)) {
                std::cerr << "error: unknown ordering '" << args[i] << "' (natural, rcm, md, amd)\n";
                return 2;
            }
        } else if (s == "--tol") {
            if (i + 1 >= args.size()) return usage();
            opt.tol = std::stod(args[++i]);
        } else if (s == "--static-pivot") {
            if (i + 1 >= args.size()) return usage();
            opt.static_pivot = std::stod(args[++i]);
        } else if (s == "--refine") {
            if (i + 1 >= args.size()) return usage();
            opt.refine = std::stoll(args[++i]);
        } else if (!s.empty() && s[0] == '-') return usage();
        else files.push_back(s);
    }
    if (files.empty()) return usage();

    if (files.size() == 1) {
        const FactorRun r = factor_one(files[0], opt, true);
        return r.ok ? 0 : 1;
    }

    // table mode
    std::cout << "ordering: " << to_string(opt.order) << (opt.postorder ? " + postorder" : "");
    if (opt.static_pivot > 0.0) std::cout << ", static pivot " << opt.static_pivot;
    if (opt.refine > 0) std::cout << ", refine <= " << opt.refine << " steps";
    std::cout << '\n';
    std::cout << std::left
              << std::setw(10) << "matrix" << std::right
              << std::setw(7) << "n" << std::setw(9) << "nnz(A)" << std::setw(10) << "nnz(L)"
              << std::setw(7) << "fill" << std::setw(10) << "flops"
              << std::setw(9) << "ord[ms]"
              << std::setw(9) << "sym[ms]" << std::setw(9) << "num[ms]" << std::setw(9) << "slv[ms]"
              << std::setw(6) << "neg" << std::setw(7) << "pert" << std::setw(11) << "min|d|"
              << std::setw(11) << "bwd0" << std::setw(4) << "it" << std::setw(11) << "bwd err"
              << std::setw(11) << "fwd err";
    if (opt.check_factor) std::cout << std::setw(11) << "fact err";
    std::cout << '\n';

    int failures = 0;
    for (const std::string& f : files) {
        FactorRun r;
        try {
            r = factor_one(f, opt, false);
        } catch (const std::exception& e) {
            std::cout << std::left << std::setw(10) << base_name(f) << "  error: " << e.what() << '\n';
            ++failures;
            continue;
        }
        std::cout << std::left << std::setw(10) << r.name << std::right
                  << std::setw(7) << r.n << std::setw(9) << r.nnz_a << std::setw(10) << r.nnz_l
                  << std::fixed << std::setprecision(2) << std::setw(7) << r.fill
                  << std::scientific << std::setprecision(2) << std::setw(10) << r.flops
                  << std::fixed << std::setprecision(1)
                  << std::setw(9) << r.t_order * 1e3
                  << std::setw(9) << r.t_sym * 1e3 << std::setw(9) << r.t_num * 1e3
                  << std::setw(9) << r.t_solve * 1e3
                  << std::setw(6) << r.neg << std::setw(7) << r.perturbed
                  << std::scientific << std::setprecision(2) << std::setw(11) << r.dmin;
        if (r.ok) {
            std::cout << std::setw(11) << r.backward0 << std::setw(4) << r.refine_iters
                      << std::setw(11) << r.backward << std::setw(11) << r.forward;
            if (opt.check_factor) {
                if (r.recon >= 0.0) std::cout << std::setw(11) << r.recon;
                else std::cout << std::setw(11) << "n/a";
            }
        } else {
            std::cout << "  FAILED: zero pivot at k = " << r.failed_pivot;
            ++failures;
        }
        std::cout << std::defaultfloat << '\n';
    }
    return failures == 0 ? 0 : 1;
}

// -- orderings: compare fill for every ordering --

int cmd_orderings(const std::vector<std::string>& args) {
    if (args.empty()) return usage();
    const Ordering methods[] = {Ordering::Natural, Ordering::RCM, Ordering::MinimumDegree, Ordering::AMD};

    std::cout << "nnz(L) below the diagonal and ordering time [ms] per method; flops of the numeric phase\n";
    std::cout << std::left << std::setw(10) << "matrix" << std::right
              << std::setw(7) << "n" << std::setw(9) << "nnz(A)";
    for (Ordering m : methods) {
        std::cout << " |" << std::setw(9) << to_string(m) << std::setw(9) << "flops" << std::setw(8) << "t[ms]";
    }
    std::cout << " |" << std::setw(9) << "amd/nat" << '\n';

    int failures = 0;
    for (const std::string& f : args) {
        CscMatrix a;
        try {
            a = read_matrix_market(f);
            if (!a.is_square()) throw InvalidMatrix("matrix is not square");
        } catch (const std::exception& e) {
            std::cout << std::left << std::setw(10) << base_name(f) << "  error: " << e.what() << '\n';
            ++failures;
            continue;
        }
        std::cout << std::left << std::setw(10) << base_name(f) << std::right
                  << std::setw(7) << a.nrows() << std::setw(9) << a.nnz();
        Index natural_nnz = 0, amd_nnz = 0;
        for (Ordering m : methods) {
            const auto t0 = Clock::now();
            const std::vector<Index> p = compute_ordering(a, m);
            const double t = seconds_since(t0);
            const FillStats fs = fill_for_ordering(a, p);
            if (m == Ordering::Natural) natural_nnz = fs.nnz_L;
            if (m == Ordering::AMD) amd_nnz = fs.nnz_L;
            std::cout << " |" << std::setw(9) << fs.nnz_L
                      << std::scientific << std::setprecision(2) << std::setw(9) << fs.flops
                      << std::fixed << std::setprecision(1) << std::setw(8) << t * 1e3 << std::defaultfloat;
        }
        std::cout << " |" << std::fixed << std::setprecision(3) << std::setw(9)
                  << (natural_nnz > 0 ? static_cast<double>(amd_nnz) / static_cast<double>(natural_nnz) : 0.0)
                  << std::defaultfloat << '\n';
    }
    return failures == 0 ? 0 : 1;
}

}  // namespace

int run_batch(const std::vector<std::string>& words) {
    if (words.empty()) return usage();
    const std::string& cmd = words[0];
    const std::vector<std::string> args(words.begin() + 1, words.end());

    if (cmd == "help" || cmd == "--help" || cmd == "-h") {
        print_usage(std::cout);
        return 0;
    }
    try {
        if (cmd == "info") return cmd_info(args);
        if (cmd == "print") return cmd_print(args);
        if (cmd == "convert") return cmd_convert(args);
        if (cmd == "factor") return cmd_factor(args);
        if (cmd == "orderings") return cmd_orderings(args);
        std::cerr << "unknown command '" << cmd << "'\n";
        return usage();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
