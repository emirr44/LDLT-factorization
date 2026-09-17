// batch.hpp - the non-interactive commands (info, print, convert, factor,
// orderings), for scripts and for measuring many matrices at once.
#pragma once

#include <string>
#include <vector>

// `words` is the command line without the program name and without the
// global options, e.g. {"factor", "a.mtx", "--check-factor"}.
// Returns the process exit code.
int run_batch(const std::vector<std::string>& words);
