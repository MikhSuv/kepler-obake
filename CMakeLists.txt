cmake_minimum_required(VERSION 3.20)

project(kepler-obake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Экспортировать compile_commands.json (для clangd).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Флаги компиляции для Debug и Release (GCC/Clang).
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

# obake-config.cmake находит obake::obake, что тянет за собой
# mp++, Boost, fmt, TBB, libbacktrace (см. медиаторы obake).
find_package(obake REQUIRED)

# Статическая библиотека: ряды z1/z2 и численная оценка.
add_library(kepler STATIC
    src/series.cpp
    src/evaluate.cpp
    src/exactkepler.cpp)
target_include_directories(kepler PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(kepler PUBLIC obake::obake)

# Демонстрационный исполняемый файл.
add_executable(kepler_z1 src/main.cpp)
target_link_libraries(kepler_z1 PRIVATE kepler)
