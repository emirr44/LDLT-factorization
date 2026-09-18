set(LDLT_NAME ldlt)

# --- SQLite, from the bundled public-domain amalgamation --------------------
add_library(sqlite3 STATIC ${CMAKE_CURRENT_LIST_DIR}/third_party/sqlite/sqlite3.c)
target_include_directories(sqlite3 PUBLIC ${CMAKE_CURRENT_LIST_DIR}/third_party/sqlite)
target_compile_definitions(sqlite3 PRIVATE SQLITE_THREADSAFE=0 SQLITE_OMIT_LOAD_EXTENSION)
if(NOT MSVC)
  target_compile_options(sqlite3 PRIVATE -w)      # third-party code: no warnings
endif()

# --- the factorization library (src/core) ----------------------------------
file(GLOB LDLT_CORE_SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/core/*.cpp)
file(GLOB LDLT_CORE_HEADERS ${CMAKE_CURRENT_LIST_DIR}/src/core/*.hpp)

add_library(ldlt_core ${LDLT_CORE_SOURCES} ${LDLT_CORE_HEADERS})
target_include_directories(ldlt_core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src/core)
if(MSVC)
  target_compile_options(ldlt_core PRIVATE /W4)
else()
  target_compile_options(ldlt_core PRIVATE -Wall -Wextra -Wpedantic -Wshadow -Wconversion)
endif()

source_group("core" FILES ${LDLT_CORE_SOURCES} ${LDLT_CORE_HEADERS})

# --- the ldlt command line tool (src/cli) ----------------------------------
# Record where the project lives, so the tool finds data/ from any directory.
set(LDLT_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(WRITE ${LDLT_GENERATED_DIR}/ldlt_config.hpp
  "#pragma once\n#define LDLT_PROJECT_ROOT \"${CMAKE_CURRENT_LIST_DIR}\"\n")

file(GLOB LDLT_CLI_SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/cli/*.cpp)
file(GLOB LDLT_CLI_HEADERS ${CMAKE_CURRENT_LIST_DIR}/src/cli/*.hpp)

add_executable(${LDLT_NAME} ${LDLT_CLI_SOURCES} ${LDLT_CLI_HEADERS})
target_include_directories(${LDLT_NAME} PRIVATE ${LDLT_GENERATED_DIR})
target_link_libraries(${LDLT_NAME} PRIVATE ldlt_core sqlite3)
if(MSVC)
  target_compile_options(${LDLT_NAME} PRIVATE /W4)
  # getenv() is standard C++; MSVC's "secure CRT" deprecation is not relevant
  target_compile_definitions(${LDLT_NAME} PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(${LDLT_NAME} PRIVATE -Wall -Wextra)
endif()

source_group("cli" FILES ${LDLT_CLI_SOURCES} ${LDLT_CLI_HEADERS})

# Opening the generated project in an IDE should let Run start ldlt directly:
# Visual Studio otherwise picks ALL_BUILD, and Xcode gets an explicit scheme
# that starts in the project folder.
set_property(DIRECTORY ${CMAKE_CURRENT_LIST_DIR} PROPERTY VS_STARTUP_PROJECT ${LDLT_NAME})
set_target_properties(${LDLT_NAME} PROPERTIES
  VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
  XCODE_GENERATE_SCHEME TRUE
  XCODE_SCHEME_WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
