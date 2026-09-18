#include "interactive.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "catalog.hpp"
#include "report.hpp"
#include "term.hpp"

#include "ordering.hpp"

namespace fs = std::filesystem;

namespace {

struct Session {
    std::string root;
    Catalog catalog;
    Ordering ordering = Ordering::AMD;
    bool has_result = false;
    FactorResult last;
};

// -- small string helpers --

std::string trim(std::string text) {
    if (text.rfind("\xEF\xBB\xBF", 0) == 0) text.erase(0, 3);
    const auto b = text.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = text.find_last_not_of(" \t\r\n");
    return text.substr(b, e - b + 1);
}

std::string lower(std::string text) {
    for (char& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

void split_command(const std::string& line, std::string& cmd, std::string& arg) {
    const auto sp = line.find_first_of(" \t");
    if (sp == std::string::npos) {
        cmd = line;
        arg.clear();
    }
    else {
        cmd = line.substr(0, sp);
        arg = trim(line.substr(sp));
    }
}

std::string clean_path(std::string p) {
    p = trim(p);
    if (p.size() >= 2 && (p.front() == '"' || p.front() == '\'') && p.back() == p.front()) {
        return p.substr(1, p.size() - 2);
    }
#ifndef _WIN32
    std::string out;
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (p[i] == '\\' && i + 1 < p.size()) ++i;
        out += p[i];
    }
    return out;
#else
    return p;
#endif
}

std::string pad(const std::string& text, std::size_t w) {
    const std::size_t have = term::width(text);
    return have >= w ? text + " " : text + std::string(w - have, ' ');
}

std::string human_size(std::uintmax_t bytes) {
    std::ostringstream os;
    const double b = static_cast<double>(bytes);
    if (b < 1024.0) os << bytes << " B";
    else if (b < 1024.0 * 1024.0) os << std::fixed << std::setprecision(1) << b / 1024.0 << " KB";
    else os << std::fixed << std::setprecision(1) << b / (1024.0 * 1024.0) << " MB";
    return os.str();
}

std::string display_path(const fs::path& p) {
    std::error_code ec;
    const fs::path rel = fs::relative(p, fs::current_path(ec), ec);
    if (!ec && !rel.empty() && rel.generic_string().rfind("..", 0) != 0) {
        return rel.generic_string();
    }
    return p.generic_string();
}

std::string alias(const CatalogEntry& e) { return "mtx_" + std::to_string(e.id); }

void error(const std::string& msg) { std::cout << term::bad(msg) << '\n'; }
void info(const std::string& msg)  { std::cout << term::dim(msg) << '\n'; }

// -- screens --

void welcome(Session& session, const std::string& notice) {
    const auto& entries = session.catalog.entries();
    const std::string dot = term::sym(" · ", " | ");

    std::vector<std::string> lines;
    lines.push_back("");
    lines.push_back(term::bold("Welcome to ldlt") + term::dim(" - sparse LDL") +
                    term::dim(term::sym("ᵀ", "^T")) + term::dim(" factorization"));
    lines.push_back("");
    lines.push_back(std::to_string(entries.size()) + " matrices in the catalog" + dot +
                    "ordering: " + to_string(session.ordering));
    lines.push_back("");
    lines.push_back(term::bold("Getting started"));

    auto item = [&lines](const std::string& cmd, const std::string& what) {
        lines.push_back("  " + pad(term::accent(cmd), 21) + term::dim(what));
    };

    if (entries.empty()) {
        item("import <file.mtx>", "add a matrix to begin");
    }
    else {
        const std::string range =
            entries.size() == 1
                ? alias(entries.front())
                : alias(entries.front()) + term::sym(" … ", " ... ") + alias(entries.back());
        item(range, "factorize a bundled matrix, smallest first");
        item("list", "show every matrix in the catalog");
        item("import <file.mtx>", "add your own matrix, or a whole folder");
    }
    item("order <method>", "natural, rcm, md or amd");
    item("save", "write the full L, D and P to a text file");
    item("help", "all commands");
    lines.push_back("");

    std::cout << '\n';
    term::box("ldlt", lines);

    if (!notice.empty()) {
        std::istringstream in(notice);
        std::string l;
        while (std::getline(in, l)) std::cout << term::warn("  " + l) << '\n';
    }
}

void help() {
    auto row = [](const std::string& cmd, const std::string& what) {
        std::cout << "  " << pad(term::accent(cmd), 24) << what << '\n';
    };
    std::cout << '\n' << term::bold("Commands") << '\n';
    row("mtx_<n>  or  <name>", "factorize a matrix from the catalog");
    row("<path/to/file.mtx>", "factorize a file directly, without importing it");
    row("list", "show the catalog");
    row("import <file|folder>", "add a .mtx file, or every .mtx in a folder");
    row("remove <mtx_n|name>", "take a matrix out of the catalog (the file is kept)");
    row("order [method]", "show or set the ordering: natural, rcm, md, amd");
    row("save [file.txt]", "write the last result in full: summary, L, D and P");
    row("clear", "clear the screen");
    row("help", "this list");
    row("exit", "quit (also: quit, q)");
    std::cout << '\n'
              << term::dim("Non-interactive use: run  ldlt help  in a shell, e.g.") << '\n'
              << term::dim("  ldlt factor data/suitesparse/bcsstk01.mtx --check-factor") << '\n';
}

void print_summary(const FactorResult& r) {
    std::istringstream in(r.summary);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("status", 0) == 0)        std::cout << (r.ok ? term::ok(line) : term::bad(line));
        else if (line.rfind("warning", 0) == 0)  std::cout << term::warn(line);
        else                                     std::cout << line;
        std::cout << '\n';
    }
}

// -- commands --

void factorize(Session& session, const std::string& path, const std::string& label) {
    info("factorizing " + label + " with " + to_string(session.ordering) + " ...");
    std::cout << std::flush;   // show that before a large matrix takes a while

    FactorResult r = factorize_matrix(path, session.ordering);
    std::cout << '\n';
    print_summary(r);

    if (r.threw) return;

    session.last = std::move(r);
    session.has_result = true;
    std::cout << '\n'
              << term::dim("type ") << term::accent("save")
              << term::dim(" to write the full L, D and P to reports/" +
                           session.last.suggested_file_name)
              << '\n';
}

void cmd_list(Session& session) {
    const auto& entries = session.catalog.entries();
    if (entries.empty()) {
        info("the catalog is empty; use: import <file.mtx>");
        return;
    }
    std::cout << '\n'
              << term::bold("  " + pad("id", 9) + pad("name", 14) + pad("n", 10) +
                            pad("nnz", 12) + "symmetry")
              << '\n';
    for (const CatalogEntry& e : entries) {
        std::cout << "  " << term::accent(pad(alias(e), 9)) << pad(e.name, 14)
                  << pad(std::to_string(e.n), 10) << pad(std::to_string(e.nnz), 12)
                  << e.symmetry << '\n';
    }
    if (!session.catalog.persistent()) {
        std::cout << term::warn("  (temporary catalog: changes are not saved)") << '\n';
    }
}

void cmd_order(Session& session, const std::string& arg) {
    if (arg.empty()) {
        std::cout << "ordering: " << term::accent(to_string(session.ordering)) << '\n'
                  << term::dim("options: natural, rcm, md, amd") << '\n';
        return;
    }
    Ordering o;
    if (!parse_ordering(arg, o)) {
        error("unknown ordering '" + arg + "'; use natural, rcm, md or amd");
        return;
    }
    session.ordering = o;
    std::cout << "ordering set to " << term::accent(to_string(o)) << '\n';
}

void cmd_import(Session& session, const std::string& arg) {
    if (arg.empty()) {
        error("usage: import <file.mtx>  or  import <folder>");
        return;
    }
    const std::string path = clean_path(arg);
    std::error_code ec;

    if (fs::is_directory(path, ec)) {
        std::string skipped;
        const int n = session.catalog.import_folder(path, skipped);
        std::cout << term::ok("imported " + std::to_string(n) + " matrices") << '\n';
        if (!skipped.empty()) {
            std::cout << term::warn("skipped:") << '\n';
            std::istringstream in(skipped);
            std::string l;
            while (std::getline(in, l)) std::cout << term::warn("  " + l) << '\n';
        }
        return;
    }

    std::string err;
    int id = 0;
    if (!session.catalog.import_file(path, err, &id)) {
        error(err);
        return;
    }
    const CatalogEntry* e = session.catalog.find(std::to_string(id));
    std::cout << term::ok("added " + (e ? e->name : path) + " as mtx_" + std::to_string(id))
              << '\n';
    if (!session.catalog.persistent()) {
        info("(temporary catalog: it will be forgotten on exit)");
    }
}

void cmd_remove(Session& session, const std::string& arg) {
    if (arg.empty()) {
        error("usage: remove <mtx_n|name>");
        return;
    }
    const CatalogEntry* e = session.catalog.find(arg);
    if (!e) {
        error("no matrix '" + arg + "' in the catalog");
        return;
    }
    const std::string label = alias(*e) + " (" + e->name + ")";
    std::string err;
    if (!session.catalog.remove(e->id, err)) {
        error(err);
        return;
    }
    std::cout << "removed " << label << '\n';
    info("the .mtx file itself was not deleted");
}

void cmd_save(Session& session, const std::string& arg) {
    if (!session.has_result) {
        error("nothing to save yet; factorize a matrix first");
        return;
    }

    std::error_code ec;
    fs::path target;
    if (arg.empty()) {
        target = fs::path(session.root) / "reports" / session.last.suggested_file_name;
    }
    else {
        target = clean_path(arg);
        if (fs::is_directory(target, ec)) target /= session.last.suggested_file_name;
    }

    if (target.has_parent_path()) fs::create_directories(target.parent_path(), ec);

    std::ofstream f(target, std::ios::binary);
    if (!f) {
        error("cannot write " + display_path(target));
        return;
    }
    f << session.last.report;
    f.close();
    if (!f) {
        error("failed while writing " + display_path(target));
        return;
    }

    const auto size = fs::file_size(target, ec);
    std::cout << term::ok("saved ") << display_path(target)
              << term::dim(ec ? std::string() : "  (" + human_size(size) + ")") << '\n';
}

}  // namespace

int run_interactive(const std::string& project_root) {
    Session session;
    session.root = project_root;

    std::string notice;
    session.catalog.open(project_root + "/data/catalog.db",
                         project_root + "/data/suitesparse", project_root, notice);
    welcome(session, notice);

    std::string line;
    for (;;) {
        std::cout << '\n'
                  << term::accent("ldlt")
                  << term::dim(std::string(" (") + to_string(session.ordering) + ")")
                  << term::accent(term::sym(" › ", " > ")) << std::flush;

        if (!std::getline(std::cin, line)) {   // end of input
            std::cout << '\n';
            break;
        }
        line = trim(line);
        if (line.empty()) continue;

        std::string cmd, arg;
        split_command(line, cmd, arg);
        const std::string c = lower(cmd);

        if (c == "exit" || c == "quit" || c == "q") break;
        if (c == "help" || c == "?")         { help(); continue; }
        if (c == "list" || c == "ls")        { cmd_list(session); continue; }
        if (c == "order")                    { cmd_order(session, arg); continue; }
        if (c == "import")                   { cmd_import(session, arg); continue; }
        if (c == "remove" || c == "rm")      { cmd_remove(session, arg); continue; }
        if (c == "save")                     { cmd_save(session, arg); continue; }
        if (c == "clear" || c == "cls")      { term::clear(); welcome(session, ""); continue; }

        // not a command: a matrix from the catalog, or a .mtx file on disk
        if (const CatalogEntry* e = session.catalog.find(line)) {
            factorize(session, e->path, alias(*e) + " (" + e->name + ")");
            continue;
        }
        const std::string path = clean_path(line);
        std::error_code ec;
        if (fs::path(path).extension() == ".mtx" && fs::is_regular_file(path, ec)) {
            factorize(session, path, fs::path(path).filename().string());
            continue;
        }

        error("'" + line + "' is not a command or a matrix in the catalog");
        info("type help for commands, or list for the matrices");
    }

    std::cout << term::dim("bye") << '\n';
    return 0;
}
