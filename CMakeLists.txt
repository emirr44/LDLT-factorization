cmake_minimum_required(VERSION 3.15)
project(ldlt LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- SQLite, from the bundled public-domain amalgamation --------------------
add_library(sqlite3 STATIC third_party/sqlite/sqlite3.c)
target_include_directories(sqlite3 PUBLIC third_party/sqlite)
target_compile_definitions(sqlite3 PRIVATE SQLITE_THREADSAFE=0 SQLITE_OMIT_LOAD_EXTENSION)
if(NOT MSVC)
  target_compile_options(sqlite3 PRIVATE -w)      # third-party code: no warnings
endif()

# --- the factorization library ---------------------------------------------
add_library(ldlt_core
  src/coo.cpp
  src/csc.cpp
  src/matrix_market.cpp
  src/symbolic.cpp
  src/numeric.cpp
  src/validate.cpp
  src/etree.cpp
  src/ordering.cpp
  src/refine.cpp
)
target_include_directories(ldlt_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(MSVC)
  target_compile_options(ldlt_core PRIVATE /W4)
else()
  target_compile_options(ldlt_core PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion)
endif()

# --- the ldlt command line tool --------------------------------------------
# Record where the project lives, so the tool finds data/ from any directory.
set(LDLT_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(WRITE ${LDLT_GENERATED_DIR}/ldlt_config.hpp
  "#pragma once\n#define LDLT_PROJECT_ROOT \"${CMAKE_CURRENT_SOURCE_DIR}\"\n")

file(GLOB LDLT_CLI_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/apps/cli/*.cpp)
add_executable(ldlt ${LDLT_CLI_SOURCES})
target_include_directories(ldlt PRIVATE ${LDLT_GENERATED_DIR})
target_link_libraries(ldlt PRIVATE ldlt_core sqlite3)
if(MSVC)
  target_compile_options(ldlt PRIVATE /W4)
  # getenv() is standard C++; MSVC's "secure CRT" deprecation is not relevant
  target_compile_definitions(ldlt PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(ldlt PRIVATE -Wall -Wextra)
endif()
