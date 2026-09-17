#include "ldlt/matrix_market.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>


namespace {

std::string lowercase(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

// Remove a leading UTF-8 byte order mark.  Windows editors (Notepad, and
// PowerShell's Set-Content -Encoding utf8) prepend one; without this the mark
// becomes part of the first token and the banner check fails on a valid file.
void strip_bom(std::string& s) {
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF) {
        s.erase(0, 3);
    }
}

// Read the next line that is not blank.  Returns false at end of stream.
bool next_nonblank_line(std::istream& in, std::string& line, Index& lineno) {
    while (std::getline(in, line)) {
        ++lineno;
        if (!trim(line).empty()) return true;
    }
    return false;
}

[[noreturn]] void fail(const std::string& what, Index lineno) {
    throw IoError("Matrix Market: " + what + " (line " + std::to_string(lineno) + ")");
}

}  // namespace

std::string to_string(MatrixMarketInfo::Field f) {
    switch (f) {
        case MatrixMarketInfo::Field::Double: return "real";
        case MatrixMarketInfo::Field::Integer: return "integer";
        case MatrixMarketInfo::Field::Pattern: return "pattern";
    }
    return "?";
}

std::string to_string(MatrixMarketInfo::Symmetry s) {
    switch (s) {
        case MatrixMarketInfo::Symmetry::General: return "general";
        case MatrixMarketInfo::Symmetry::Symmetric: return "symmetric";
        case MatrixMarketInfo::Symmetry::SkewSymmetric: return "skew-symmetric";
    }
    return "?";
}

// ---------------------------------------------------------------------------
// reader
// ---------------------------------------------------------------------------

CooMatrix read_matrix_market_coo(std::istream& in, bool expand_symmetric,
                                 MatrixMarketInfo* info_out) {
    MatrixMarketInfo info;
    Index lineno = 0;
    std::string line;

    // --- header -------------------------------------------------------------
    if (!std::getline(in, line)) fail("empty input", lineno);
    ++lineno;
    strip_bom(line);
    {
        std::istringstream hs(lowercase(trim(line)));
        std::string banner, object, format, field, symmetry;
        hs >> banner >> object >> format >> field >> symmetry;
        if (banner != "%%matrixmarket") fail("missing %%MatrixMarket banner", lineno);
        if (object != "matrix") fail("unsupported object '" + object + "'", lineno);
        if (format != "coordinate") {
            fail("unsupported format '" + format + "' (only coordinate is supported)", lineno);
        }
        if (field == "real") info.field = MatrixMarketInfo::Field::Double;
        else if (field == "integer") info.field = MatrixMarketInfo::Field::Integer;
        else if (field == "pattern") info.field = MatrixMarketInfo::Field::Pattern;
        else fail("unsupported field '" + field + "' (real, integer, pattern supported)", lineno);

        if (symmetry == "general") info.symmetry = MatrixMarketInfo::Symmetry::General;
        else if (symmetry == "symmetric") info.symmetry = MatrixMarketInfo::Symmetry::Symmetric;
        else if (symmetry == "skew-symmetric") info.symmetry = MatrixMarketInfo::Symmetry::SkewSymmetric;
        else fail("unsupported symmetry '" + symmetry + "'", lineno);
    }

    // --- comments and size line ----------------------------------------------
    for (;;) {
        if (!next_nonblank_line(in, line, lineno)) fail("missing size line", lineno);
        const std::string t = trim(line);
        if (t[0] == '%') {
            info.comments.push_back(t.substr(1));
            continue;
        }
        std::istringstream ss(t);
        long long m = 0, n = 0, nz = 0;
        if (!(ss >> m >> n >> nz)) fail("malformed size line '" + t + "'", lineno);
        std::string extra;
        if (ss >> extra) fail("unexpected token on size line: '" + extra + "'", lineno);
        if (m < 0 || n < 0 || nz < 0) fail("negative size", lineno);
        info.nrows = static_cast<Index>(m);
        info.ncols = static_cast<Index>(n);
        info.nnz_in_file = static_cast<Index>(nz);
        break;
    }
    if (info.is_symmetric_storage() && info.nrows != info.ncols) {
        fail("symmetric matrix must be square", lineno);
    }

    // --- entries --------------------------------------------------------------
    const bool skew = info.symmetry == MatrixMarketInfo::Symmetry::SkewSymmetric;
    const bool expand = expand_symmetric && info.is_symmetric_storage();

    CooMatrix coo(info.nrows, info.ncols);
    coo.reserve(expand ? 2 * info.nnz_in_file : info.nnz_in_file);

    for (Index k = 0; k < info.nnz_in_file; ++k) {
        if (!next_nonblank_line(in, line, lineno)) {
            fail("file ends after " + std::to_string(k) + " of " +
                     std::to_string(info.nnz_in_file) + " entries", lineno);
        }
        std::istringstream ss(line);
        long long i = 0, j = 0;
        if (!(ss >> i >> j)) fail("malformed entry '" + trim(line) + "'", lineno);

        Real v = 1.0;
        if (info.field != MatrixMarketInfo::Field::Pattern) {
            if (!(ss >> v)) fail("missing value in entry '" + trim(line) + "'", lineno);
        }
        if (i < 1 || i > info.nrows || j < 1 || j > info.ncols) {
            fail("index (" + std::to_string(i) + ", " + std::to_string(j) + ") out of range",
                 lineno);
        }
        const Index r = static_cast<Index>(i - 1);
        const Index c = static_cast<Index>(j - 1);
        if (info.is_symmetric_storage() && r < c) {
            fail("symmetric file must store the lower triangle only, found (" +
                     std::to_string(i) + ", " + std::to_string(j) + ")", lineno);
        }
        if (skew && r == c) {
            fail("skew-symmetric matrix cannot have diagonal entries", lineno);
        }
        coo.add(r, c, v);
        if (expand && r != c) coo.add(c, r, skew ? -v : v);
    }

    if (info_out) *info_out = info;
    return coo;
}

CooMatrix read_matrix_market_coo(const std::string& path, bool expand_symmetric,
                                 MatrixMarketInfo* info) {
    std::ifstream f(path);
    if (!f) throw IoError("Matrix Market: cannot open '" + path + "'");
    try {
        return read_matrix_market_coo(f, expand_symmetric, info);
    } catch (const IoError& e) {
        throw IoError(std::string(e.what()) + " in '" + path + "'");
    }
}

CscMatrix read_matrix_market(std::istream& in, bool expand_symmetric, MatrixMarketInfo* info) {
    return CscMatrix::from_coo(read_matrix_market_coo(in, expand_symmetric, info));
}

CscMatrix read_matrix_market(const std::string& path, bool expand_symmetric,
                             MatrixMarketInfo* info) {
    return CscMatrix::from_coo(read_matrix_market_coo(path, expand_symmetric, info));
}

// ---------------------------------------------------------------------------
// writer
// ---------------------------------------------------------------------------

void write_matrix_market(std::ostream& out, const CscMatrix& a, bool symmetric_storage,
                         const std::string& comment) {
    if (symmetric_storage && !a.is_square()) {
        throw InvalidMatrix("write_matrix_market: symmetric storage needs a square matrix");
    }
    CscMatrix lower;
    if (symmetric_storage) lower = a.lower_triangle(true);
    const CscMatrix& m = symmetric_storage ? lower : a;

    out << "%%MatrixMarket matrix coordinate real "
        << (symmetric_storage ? "symmetric" : "general") << '\n';
    if (!comment.empty()) {
        std::istringstream cs(comment);
        std::string cl;
        while (std::getline(cs, cl)) out << '%' << cl << '\n';
    }
    out << m.nrows() << ' ' << m.ncols() << ' ' << m.nnz() << '\n';

    std::ios old_state(nullptr);
    old_state.copyfmt(out);
    out << std::setprecision(17);
    const auto& colptr = m.colptr();
    const auto& rowind = m.rowind();
    const auto& values = m.values();
    for (Index j = 0; j < m.ncols(); ++j) {
        for (Index p = colptr[static_cast<std::size_t>(j)];
             p < colptr[static_cast<std::size_t>(j) + 1]; ++p) {
            out << rowind[static_cast<std::size_t>(p)] + 1 << ' ' << j + 1 << ' '
                << values[static_cast<std::size_t>(p)] << '\n';
        }
    }
    out.copyfmt(old_state);
    if (!out) throw IoError("Matrix Market: write failed");
}

void write_matrix_market(const std::string& path, const CscMatrix& a, bool symmetric_storage,
                         const std::string& comment) {
    std::ofstream f(path);
    if (!f) throw IoError("Matrix Market: cannot create '" + path + "'");
    write_matrix_market(f, a, symmetric_storage, comment);
}

