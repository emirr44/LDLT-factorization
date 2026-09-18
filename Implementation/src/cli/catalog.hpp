// catalog.hpp - the SQLite database of matrices the tool knows about.

#pragma once

#include <string>
#include <vector>

struct sqlite3;


struct CatalogEntry {
    int id = 0;
    std::string name;        // file name without .mtx
    std::string path;
    long long n = 0;
    long long nnz = 0;
    std::string field;       // real | integer | pattern
    std::string symmetry;    // general | symmetric | skew-symmetric
};


class Catalog {
public:
    Catalog() = default;
    ~Catalog();
    Catalog(const Catalog&) = delete;
    Catalog& operator=(const Catalog&) = delete;

    void open(const std::string& db_file, const std::string& seed_dir,
              const std::string& project_root, std::string& notice);

    const std::vector<CatalogEntry>& entries() const { return entries_; }

    bool persistent() const { return persistent_; }

    // accepts "mtx_3", "mtx3", "3", or a matrix name
    const CatalogEntry* find(const std::string& key) const;

    bool import_file(const std::string& path, std::string& err, int* out_id = nullptr);
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
