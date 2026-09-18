# Sparse LDLᵀ factorization

C++ implementation of the sparse LDL' factorization for symmetric positive
definite and indefinite matrices, following Timothy Davis, *Direct Methods for
Sparse Linear Systems*, and the educational `LDL` package in SuiteSparse.

It comes with `ldlt`, a command-line tool backed by a small SQLite catalog of
test matrices.

## Quick start

Download or clone repository, navigate to the Implementation file in Terminal with Admin privileges (cd path_to_project_src), and enter following commands to start CLI tool

**macOS / Linux**

```bash
bash build.sh (or powershell -ExecutionPolicy Bypass -File .\build.ps1)
./build/ldlt
```

**Windows**

```powershell
.\build.ps1
.\build\ldlt.exe
```

Alternatively, open it with CMake app, select ldlt as Startup Project and compile directly. That will open an app.

The first build also compiles SQLite, which takes about half a minute; later
builds reuse it. Nothing has to be installed because SQLite is bundled as source.

## Using ldlt

Started without arguments, `ldlt` opens an interactive session:

```
╭─ ldlt ──────────────────────────────────────────────────────────────╮
│                                                                     │
│  Welcome to ldlt - sparse LDLᵀ factorization                        │
│                                                                     │
│  18 matrices in the catalog · ordering: amd                         │
│                                                                     │
│  Getting started                                                    │
│    mtx_1 … mtx_18       factorize a bundled matrix, smallest first  │
│    list                 show every matrix in the catalog            │
│    import <file.mtx>    add your own matrix, or a whole folder      │
│    order <method>       natural, rcm, md or amd                     │
│    save                 write the full L, D and P to a text file    │
│    help                 all commands                                │
│                                                                     │
╰─────────────────────────────────────────────────────────────────────╯

ldlt (amd) › mtx_1
```

| Command | Effect |
|---|---|
| `mtx_<n>` or a name | factorize a matrix from the catalog (`mtx_3`, `bcsstk03`) |
| `<path/to/file.mtx>` | factorize a file directly, without importing it |
| `list` | show the catalog |
| `import <file\|folder>` | add a `.mtx` file, or every `.mtx` in a folder |
| `remove <mtx_n\|name>` | take a matrix out of the catalog (the file is kept) |
| `order [method]` | show or set the fill-reducing ordering |
| `save [file.txt]` | write the last result in full — summary, L, D and P — by default to `reports/` |
| `help`, `clear`, `exit` | |

Each factorization prints n, nnz(A), nnz(L), fill ratio, flop count, the time
of every phase, the inertia (pivot signs), static-pivoting statistics, and the
backward and forward error before and after iterative refinement.

Global options: `--plain` turns off colour and box drawing, `--root <dir>`
points at a project folder containing `data/`.

### The matrix catalog

The catalog is an SQLite database, `data/catalog.db`. It is generated, not
committed. Paths are stored relative to the project, so
the project can be moved or copied freely. If the database file cannot be
created, the tool keeps working with a temporary in-memory catalog.


## Test matrices

The 18 matrices used throughout are committed in `data/suitesparse/`, about
7 MB, so the project is self-contained. They come from the
[SuiteSparse Matrix Collection](https://sparse.tamu.edu/):

- `bcsstk01`–`bcsstk06`, `bcsstk08`–`bcsstk11`, `bcsstk13`, `bcsstk14` from the
  **HB** (Harwell-Boeing) group: symmetric positive definite stiffness matrices
- `laser`, `qpband`, `aug3d`, `tuma2`, `sit100`, `stokes64` from the
  **GHS_indef** group: symmetric indefinite systems

Any other Matrix Market file can be added with `import` in the interactive
session, or by copying it into `data/suitesparse/`.


/dense kernels, dense-row handling in AMD.
