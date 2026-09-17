#include "catalog.hpp"

#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <sstream>

#include "sqlite3.h"

#include "ldlt/matrix_market.hpp"

namespace fs = std::filesystem;

namespace {

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string base_name(const std::string& path) {
    std::string b = fs::path(path).filename().string();
    const auto d = b.rfind(".mtx");
    if (d != std::string::npos && d + 4 == b.size()) b.erase(d);
    return b;
}

std::string column_text(sqlite3_stmt* st, int col) {
    const unsigned char* t = sqlite3_column_text(st, col);
    return t ? std::string(reinterpret_cast<const char*>(t)) : std::string();
}

// Reads a Matrix Market file completely. That is the validation: a matrix the
// tool could not factorize never enters the catalog.
bool inspect(const std::string& path, CatalogEntry& out, std::string& err) {
    try {
        MatrixMarketInfo info;
        const CscMatrix a = read_matrix_market(path, true, &info);
        if (!a.is_square()) {
            err = base_name(path) + ": matrix is not square";
            return false;
        }
        out.name = base_name(path);
        out.path = path;
        out.n = a.nrows();
        out.nnz = a.nnz();
        out.field = to_string(info.field);
        out.symmetry = to_string(info.symmetry);
        return true;
    }
    catch (const std::exception& e) {
        err = base_name(path) + ": " + e.what();
        return false;
    }
}

std::vector<std::string> mtx_files_in(const std::string& folder) {
    std::vector<std::string> files;
    std::error_code ec;
    if (!fs::is_directory(folder, ec)) return files;
    for (const auto& de : fs::directory_iterator(folder, ec)) {
        if (de.is_regular_file(ec) && de.path().extension() == ".mtx") {
            files.push_back(de.path().generic_string());
        }
    }
    return files;
}

// Runs one statement that returns no rows.
bool exec(sqlite3* db, const char* sql, std::string& err) {
    char* msg = nullptr;
    if (sqlite3_exec(db, sql, nullptr, nullptr, &msg) != SQLITE_OK) {
        err = msg ? msg : sqlite3_errmsg(db);
        sqlite3_free(msg);
        return false;
    }
    return true;
}

}  // namespace


Catalog::~Catalog() {
    if (db_) sqlite3_close(db_);
}


// ---------------------------------------------------------------------------
// opening
// ---------------------------------------------------------------------------

bool Catalog::open_database(const std::string& db_file, std::string& err) {
    std::error_code ec;
    fs::create_directories(fs::path(db_file).parent_path(), ec);
    if (ec) {
        err = ec.message();
        return false;
    }
    if (sqlite3_open_v2(db_file.c_str(), &db_,
                        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, nullptr) != SQLITE_OK) {
        err = db_ ? sqlite3_errmsg(db_) : "out of memory";
        sqlite3_close(db_);
        db_ = nullptr;
        return false;
    }
    return true;
}


bool Catalog::create_schema(std::string& err) {
    return exec(db_,
                "CREATE TABLE IF NOT EXISTS Matrix ("
                "  ID    INTEGER PRIMARY KEY AUTOINCREMENT,"
                "  Name  TEXT    NOT NULL,"
                "  Path  TEXT    NOT NULL UNIQUE,"
                "  N     INTEGER NOT NULL,"
                "  NNZ   INTEGER NOT NULL,"
                "  Field TEXT    NOT NULL,"
                "  Symm  TEXT    NOT NULL)",
                err);
}


void Catalog::open(const std::string& db_file, const std::string& seed_dir,
                   const std::string& project_root, std::string& notice) {
    root_ = fs::path(project_root).generic_string();
    notice.clear();

    std::string err;
    if (!open_database(db_file, err) || !create_schema(err)) {
        // Fall back to a database that lives only in memory. Every command
        // still works; imports are simply forgotten on exit.
        if (db_) sqlite3_close(db_);
        db_ = nullptr;
        persistent_ = false;
        sqlite3_open(":memory:", &db_);
        std::string ignored;
        create_schema(ignored);
        notice = "could not open " + db_file + " (" + err + ")\n"
                 "using a temporary catalog; imports will not be kept";
    }

    reload(err);

    if (entries_.empty()) {
        std::string skipped;
        if (import_folder(seed_dir, skipped) == 0) {
            if (!notice.empty()) notice += "\n";
            notice += "no matrices found in " + seed_dir + "; use: import <file.mtx>";
        }
    }
}


// ---------------------------------------------------------------------------
// paths: stored relative to the project root when the file is inside it
// ---------------------------------------------------------------------------

std::string Catalog::to_stored_path(const std::string& absolute) const {
    std::error_code ec;
    const fs::path rel = fs::relative(absolute, root_, ec);
    if (ec || rel.empty()) return absolute;
    const std::string s = rel.generic_string();
    if (s.rfind("..", 0) == 0) return absolute;   // outside the project
    return s;
}


std::string Catalog::from_stored_path(const std::string& stored) const {
    const fs::path p(stored);
    if (p.is_absolute()) return p.generic_string();
    return (fs::path(root_) / p).generic_string();
}


// ---------------------------------------------------------------------------
// reading
// ---------------------------------------------------------------------------

bool Catalog::reload(std::string& err) {
    entries_.clear();

    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db_,
                           "SELECT ID, Name, Path, N, NNZ, Field, Symm FROM Matrix ORDER BY ID",
                           -1, &st, nullptr) != SQLITE_OK) {
        err = sqlite3_errmsg(db_);
        return false;
    }

    while (sqlite3_step(st) == SQLITE_ROW) {
        CatalogEntry e;
        e.id = sqlite3_column_int(st, 0);
        e.name = column_text(st, 1);
        e.path = from_stored_path(column_text(st, 2));
        e.n = sqlite3_column_int64(st, 3);
        e.nnz = sqlite3_column_int64(st, 4);
        e.field = column_text(st, 5);
        e.symmetry = column_text(st, 6);

        // A row whose file is gone would only fail later, so leave it out.
        std::error_code ec;
        if (!fs::exists(e.path, ec)) continue;

        entries_.push_back(std::move(e));
    }
    sqlite3_finalize(st);
    return true;
}


const CatalogEntry* Catalog::find(const std::string& raw) const {
    const std::string key = lower(raw);

    // "mtx_3", "mtx3" and "3" all mean id 3
    std::string digits = key;
    if (key.rfind("mtx_", 0) == 0)     digits = key.substr(4);
    else if (key.rfind("mtx", 0) == 0) digits = key.substr(3);

    const bool numeric = !digits.empty() && digits.size() <= 9 &&
                         std::all_of(digits.begin(), digits.end(),
                                     [](char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; });
    if (numeric) {
        const int id = std::stoi(digits);
        for (const CatalogEntry& e : entries_) {
            if (e.id == id) return &e;
        }
        return nullptr;
    }

    for (const CatalogEntry& e : entries_) {
        if (lower(e.name) == key) return &e;
    }
    return nullptr;
}


// ---------------------------------------------------------------------------
// writing
// ---------------------------------------------------------------------------

bool Catalog::insert(const CatalogEntry& e, std::string& err, int* out_id) {
    const std::string stored = to_stored_path(e.path);

    // Insert, or refresh the metadata if this path is already listed. The
    // upsert keeps the existing id, so mtx_<id> never changes under the user.
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db_,
                           "INSERT INTO Matrix (Name, Path, N, NNZ, Field, Symm)"
                           " VALUES (?1, ?2, ?3, ?4, ?5, ?6)"
                           " ON CONFLICT(Path) DO UPDATE SET"
                           "   Name = excluded.Name, N = excluded.N, NNZ = excluded.NNZ,"
                           "   Field = excluded.Field, Symm = excluded.Symm",
                           -1, &st, nullptr) != SQLITE_OK) {
        err = sqlite3_errmsg(db_);
        return false;
    }
    // SQLITE_STATIC: the strings outlive the statement, so no copy is needed.
    sqlite3_bind_text(st, 1, e.name.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_text(st, 2, stored.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int64(st, 3, e.n);
    sqlite3_bind_int64(st, 4, e.nnz);
    sqlite3_bind_text(st, 5, e.field.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_text(st, 6, e.symmetry.c_str(), -1, SQLITE_STATIC);

    const bool done = sqlite3_step(st) == SQLITE_DONE;
    if (!done) err = sqlite3_errmsg(db_);
    sqlite3_finalize(st);
    if (!done) return false;

    if (out_id) {
        *out_id = 0;
        if (sqlite3_prepare_v2(db_, "SELECT ID FROM Matrix WHERE Path = ?1",
                               -1, &st, nullptr) == SQLITE_OK) {
            sqlite3_bind_text(st, 1, stored.c_str(), -1, SQLITE_STATIC);
            if (sqlite3_step(st) == SQLITE_ROW) *out_id = sqlite3_column_int(st, 0);
            sqlite3_finalize(st);
        }
    }
    return true;
}


bool Catalog::import_file(const std::string& path, std::string& err, int* out_id) {
    std::error_code ec;
    const fs::path abs = fs::absolute(path, ec);
    if (ec || !fs::is_regular_file(abs, ec)) {
        err = path + ": no such file";
        return false;
    }

    CatalogEntry e;
    if (!inspect(abs.generic_string(), e, err)) return false;
    if (!insert(e, err, out_id)) return false;

    std::string ignored;
    reload(ignored);
    return true;
}


int Catalog::import_folder(const std::string& folder, std::string& err) {
    std::ostringstream problems;

    std::vector<CatalogEntry> found;
    for (const std::string& f : mtx_files_in(folder)) {
        CatalogEntry e;
        std::string why;
        if (inspect(f, e, why)) found.push_back(std::move(e));
        else problems << why << '\n';
    }

    // Smallest first, so the lowest ids are the quickest to factorize.
    std::sort(found.begin(), found.end(), [](const CatalogEntry& a, const CatalogEntry& b) {
        return a.n != b.n ? a.n < b.n : a.name < b.name;
    });

    int added = 0;
    std::string ignored;
    exec(db_, "BEGIN", ignored);   // one transaction: far faster than one per row
    for (const CatalogEntry& e : found) {
        std::string why;
        if (insert(e, why, nullptr)) ++added;
        else problems << e.name << ": " << why << '\n';
    }
    exec(db_, "COMMIT", ignored);

    reload(ignored);
    err = problems.str();
    return added;
}


bool Catalog::remove(int id, std::string& err) {
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(db_, "DELETE FROM Matrix WHERE ID = ?1", -1, &st, nullptr) != SQLITE_OK) {
        err = sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_int(st, 1, id);
    const bool done = sqlite3_step(st) == SQLITE_DONE;
    if (!done) err = sqlite3_errmsg(db_);
    sqlite3_finalize(st);
    if (!done) return false;

    if (sqlite3_changes(db_) == 0) {
        err = "no matrix with id " + std::to_string(id);
        return false;
    }
    std::string ignored;
    reload(ignored);
    return true;
}
