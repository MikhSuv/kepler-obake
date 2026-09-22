cmake_minimum_required(VERSION 3.20)

project(kepler-obake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Экспортировать compile_commands.json (для clangd).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# obake-config.cmake находит obake::obake, что тянет за собой
# mp++, Boost, fmt, TBB, libbacktrace (см. медиаторы obake).
find_package(obake REQUIRED)

# Статическая библиотека: ряды z1/z2 и численная оценка.
add_library(kepler STATIC
    src/series.cpp
    src/evaluate.cpp
    src/exactkepler.cpp)
target_include_directories(kepler PUBLIC src)
target_link_libraries(kepler PUBLIC obake::obake)

# Демонстрационный исполняемый файл.
add_executable(kepler_z1 src/main.cpp)
target_link_libraries(kepler_z1 PRIVATE kepler)
