
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace term {

void init(bool force_plain);

bool fancy();

const char* sym(const char* unicode, const char* ascii);

std::string accent(const std::string& s);   // orange, for the banner and prompt
std::string ok(const std::string& s);       // green
std::string bad(const std::string& s);      // red
std::string warn(const std::string& s);     // yellow
std::string dim(const std::string& s);      // grey
std::string bold(const std::string& s);

std::size_t width(const std::string& s);

void box(const std::string& title, const std::vector<std::string>& lines);

void clear();   

}
