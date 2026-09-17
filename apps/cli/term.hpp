// term.hpp - small terminal helpers: colour, UTF-8, and a framed panel.
//
// Everything degrades to plain ASCII without colour when the output is not a
// terminal (e.g. redirected into a file), when NO_COLOR is set, or when the
// user passes --plain. The program never depends on the fancy output.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace term {

// Call once at startup.
void init(bool force_plain);

// True when colour and Unicode box drawing are in use.
bool fancy();

// Pick the Unicode or the ASCII spelling of a symbol.
const char* sym(const char* unicode, const char* ascii);

// Colour wrappers; return the text unchanged when fancy() is false.
std::string accent(const std::string& s);   // orange, for the banner and prompt
std::string ok(const std::string& s);       // green
std::string bad(const std::string& s);      // red
std::string warn(const std::string& s);     // yellow
std::string dim(const std::string& s);      // grey
std::string bold(const std::string& s);

// On-screen width: counts UTF-8 code points and skips colour escapes.
std::size_t width(const std::string& s);

// Draws a framed panel with `title` in the top border.
void box(const std::string& title, const std::vector<std::string>& lines);

// Clears the screen (no-op in plain mode).
void clear();

}  // namespace term
