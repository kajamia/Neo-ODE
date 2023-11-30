cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(Neo-CLA/src)
add_subdirectory (src)

add_executable (test_ode demos/test_ode.cc)
target_link_libraries (test_ode PUBLIC ngbla)

add_executable (test_newmark demos/test_newmark.cc)
target_link_libraries (test_newmark PUBLIC ngbla)

add_executable (test_alpha demos/test_alpha.cc)
target_link_libraries (test_alpha PUBLIC ngbla)


# add_subdirectory (mass_spring)
