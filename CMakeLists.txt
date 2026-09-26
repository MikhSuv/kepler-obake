cmake_minimum_required(VERSION 3.20)

project(kepler-obake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile_commands.json (for clangd).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Debug and Release compilation flags (GCC/Clang).
if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/W4 /permissive- /Z7")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/fsanitize=address")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
else()
    set(CMAKE_CXX_FLAGS_DEBUG
        "-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wcast-align -Wformat=2 -Werror -O0 -fno-omit-frame-pointer -fsanitize=address,undefined")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address,undefined")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
endif()

# obake-config.cmake provides obake::obake, which pulls in mp++, Boost, fmt,
# TBB and libbacktrace (see obake's mediators).
find_package(obake REQUIRED)

# Static library: z1/z2 series and numerical evaluation.
add_library(kepler STATIC
    src/series.cpp
    src/evaluate.cpp
    src/exactkepler.cpp)
target_include_directories(kepler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(kepler PUBLIC obake::obake)

# Demo executable.
add_executable(kepler_z1 src/main.cpp)
target_link_libraries(kepler_z1 PRIVATE kepler)
