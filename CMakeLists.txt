cmake_minimum_required(VERSION 3.16)

project(kepler-obake LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Экспортировать compile_commands.json (для clangd).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# obake-config.cmake находит obake::obake, что тянет за собой
# mp++, Boost, fmt, TBB, libbacktrace (см. медиаторы obake).
find_package(obake REQUIRED)

add_executable(kepler_z1 src/main.cpp)
target_link_libraries(kepler_z1 PRIVATE obake::obake)