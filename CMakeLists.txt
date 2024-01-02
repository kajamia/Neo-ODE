cmake_minimum_required(VERSION 3.20)
project(Neo_ODE)

set (CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(Neo-CLA/src)
add_subdirectory (src)

add_executable (test_ode demos/test_ode.cc)
# target_link_libraries (test_ode PUBLIC ngbla)

# set(CMAKE_BUILD_TYPE Debug)


# test_exponential from test_ode for Python

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(ode demos/bind_test_ode.cc)
install (TARGETS ode DESTINATION Neoode)


add_executable (test_newmark demos/test_newmark.cc)
# target_link_libraries (test_newmark PUBLIC ngbla)

add_executable (test_alpha demos/test_alpha.cc)
# target_link_libraries (test_alpha PUBLIC ngbla)


add_subdirectory (mass_spring)
