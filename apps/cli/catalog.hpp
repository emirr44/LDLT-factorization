// catalog.hpp - the SQLite database of matrices the tool knows about.
//
// One row per Matrix Market file: an id (shown as mtx_<id>), the name, the
// path, and the metadata the list command prints. The matrix itself stays in
// its .mtx file and is re-read on every factorization.
//
//   - First run: the catalog seeds itself from data/suitesparse, smallest
//     matrix first, so mtx_1 is always a quick one.
//   - Paths are stored relative to the project root, so the database keeps
//     working when the project is moved or copied to another machine.
//   - Rows whose file has disappeared are skipped when loading.
//   - If the database file cannot be created (e.g. a read-only folder), the
//     catalog runs in memory instead: everything works, imports just are not
//     remembered after exit.
#pragma once

#include <string>
#include <vector>

struct sqlite3;


struct CatalogEntry {
    int id = 0;
    std::string name;        // file name without .mtx
    std::string path;        // absolute path, known to exist
    long long n = 0;         // rows (== columns)
    long long nnz = 0;       // stored entries after symmetric expansion
    std::string field;       // real | integer | pattern
    std::string symmetry;    // general | symmetric | skew-symmetric
};


class Catalog {
public:
    Catalog() = default;
    ~Catalog();
    Catalog(const Catalog&) = delete;
    Catalog& operator=(const Catalog&) = delete;

    // Opens (or creates) the database and loads it. Never leaves the catalog
    // unusable; `notice` explains anything the user should know about.
    void open(const std::string& db_file, const std::string& seed_dir,
              const std::string& project_root, std::string& notice);

    const std::vector<CatalogEntry>& entries() const { return entries_; }

    // False when running in memory because the file could not be opened.
    bool persistent() const { return persistent_; }

    // Accepts "mtx_3", "mtx3", "3", or a matrix name (case-insensitive).
    const CatalogEntry* find(const std::string& key) const;

    // Validates the file by reading it, then adds it (or refreshes it if it is
    // already there). `out_id` receives the id it is listed under.
    bool import_file(const std::string& path, std::string& err, int* out_id = nullptr);

    // Imports every .mtx directly inside `folder`, smallest first.
    // Returns how many were added; `err` lists the files that were rejected.
    int import_folder(const std::string& folder, std::string& err);

    bool remove(int id, std::string& err);

private:
    sqlite3* db_ = nullptr;
    bool persistent_ = true;
    std::string root_;
    std::vector<CatalogEntry> entries_;

    bool open_database(const std::string& db_file, std::string& err);
    bool create_schema(std::string& err);
    bool reload(std::string& err);
    bool insert(const CatalogEntry& e, std::string& err, int* out_id);

    std::string to_stored_path(const std::string& absolute) const;
    std::string from_stored_path(const std::string& stored) const;
};
