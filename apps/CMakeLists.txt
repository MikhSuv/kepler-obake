add_executable(z1_z2 z1_z2.cpp)
target_link_libraries(z1_z2 PRIVATE kepler)

add_executable(z3 z3.cpp)
target_link_libraries(z3 PRIVATE kepler)

add_executable(z4_z5 z4_z5.cpp)
target_link_libraries(z4_z5 PRIVATE kepler)

add_executable(tabulate evaluate.cpp)
target_link_libraries(tabulate PRIVATE kepler)
