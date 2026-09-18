#include "term.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <io.h>
#include <windows.h>
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#else
#include <unistd.h>
#endif

namespace term {

namespace {

bool g_fancy = false;

std::string wrap(const char* code, const std::string& s) {
    if (!g_fancy) return s;
    return std::string("\x1b[") + code + "m" + s + "\x1b[0m";
}

std::string repeat(const char* piece, std::size_t count) {
    std::string out;
    for (std::size_t i = 0; i < count; ++i) out += piece;
    return out;
}

bool stdout_is_terminal() {
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

}

void init(bool force_plain) {
    const bool tty = stdout_is_terminal();
    const bool forced = std::getenv("FORCE_COLOR") != nullptr;
    g_fancy = !force_plain && std::getenv("NO_COLOR") == nullptr && (tty || forced);

#ifdef _WIN32
    if (g_fancy && tty) {
        SetConsoleOutputCP(CP_UTF8);
        HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD mode = 0;
        if (h == INVALID_HANDLE_VALUE || !GetConsoleMode(h, &mode) ||
            !SetConsoleMode(h, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING)) {
            g_fancy = false;
        }
    }
#endif
}

bool fancy() { return g_fancy; }

const char* sym(const char* unicode, const char* ascii) { return g_fancy ? unicode : ascii; }

std::string accent(const std::string& s) { return wrap("38;5;209", s); }
std::string ok(const std::string& s)     { return wrap("32", s); }
std::string bad(const std::string& s)    { return wrap("31", s); }
std::string warn(const std::string& s)   { return wrap("33", s); }
std::string dim(const std::string& s)    { return wrap("2", s); }
std::string bold(const std::string& s)   { return wrap("1", s); }

std::size_t width(const std::string& s) {
    std::size_t w = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == 0x1b) {                       
            while (i < s.size() && s[i] != 'm') ++i;
            continue;
        }
        if ((c & 0xC0) != 0x80) ++w;           
    }
    return w;
}

void box(const std::string& title, const std::vector<std::string>& lines) {
    const char* h  = sym("─", "-");
    const char* v  = sym("│", "|");
    const char* tl = sym("╭", "+");
    const char* tr = sym("╮", "+");
    const char* bl = sym("╰", "+");
    const char* br = sym("╯", "+");

    const std::size_t pad = 2;
    std::size_t inner = width(title) + 4;
    for (const std::string& l : lines) inner = std::max(inner, width(l) + 2 * pad);

    std::cout << accent(std::string(tl) + h + " ") << bold(accent(title))
              << accent(" " + repeat(h, inner - width(title) - 3) + tr) << '\n';

    for (const std::string& l : lines) {
        std::cout << accent(v) << std::string(pad, ' ') << l
                  << std::string(inner - pad - width(l), ' ') << accent(v) << '\n';
    }

    std::cout << accent(std::string(bl) + repeat(h, inner) + br) << '\n';
}

void clear() {
    if (g_fancy) std::cout << "\x1b[2J\x1b[H" << std::flush;
}

}
