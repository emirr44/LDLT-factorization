// interactive.hpp - the interactive session started by running `ldlt` alone.
#pragma once

#include <string>

// Opens the matrix catalog under `project_root/data`, shows the welcome
// panel, and reads commands until exit or end of input. Returns the process
// exit code.
int run_interactive(const std::string& project_root);
