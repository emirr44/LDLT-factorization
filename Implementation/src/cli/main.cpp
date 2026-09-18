

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "batch.hpp"
#include "interactive.hpp"
#include "term.hpp"


#if __has_include("ldlt_config.hpp")
#include "ldlt_config.hpp"
#endif

namespace fs = std::filesystem;

namespace {

bool has_data(const fs::path& dir) {
    std::error_code ec;
    return fs::is_directory(dir / "data", ec);
}

std::string project_root(const std::string& from_flag) {
    if (!from_flag.empty()) return fs::path(from_flag).generic_string();
    if (const char* env = std::getenv("LDLT_ROOT")) return fs::path(env).generic_string();
#ifdef LDLT_PROJECT_ROOT
    if (has_data(LDLT_PROJECT_ROOT)) return fs::path(LDLT_PROJECT_ROOT).generic_string();
#endif
    std::error_code ec;
    return fs::current_path(ec).generic_string();
}

}  // namespace

int main(int argc, char** argv) {
    bool plain = false;
    std::string root_flag;
    std::vector<std::string> words;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--plain") {
            plain = true;
        }
        else if (a == "--root") {
            if (i + 1 >= argc) {
                std::cerr << "--root needs a folder\n";
                return 2;
            }
            root_flag = argv[++i];
        }
        else {
            words.push_back(a);
        }
    }

    term::init(plain);

    if (words.empty()) {
        const std::string root = project_root(root_flag);
        if (!has_data(root)) {
            std::cerr << "cannot find the data/ folder under " << root << "\n"
                      << "run ldlt from the project folder, or pass --root <project folder>\n";
            return 1;
        }
        return run_interactive(root);
    }
    return run_batch(words);
}
