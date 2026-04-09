cmake_minimum_required(VERSION 3.15)

project(ninteg Fortran)

set(MIN_FORTRAN_VERSION "9.3.0")

set(CMAKE_LINK_DEPENDS_USE_LINKER FALSE)

set(CMAKE_Fortran_COMPILER "gfortran")

# Set compiler flags for different build types
set(CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -fopenmp -fbacktrace")
set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -fbacktrace -ffpe-trap=overflow,zero,invalid -O0 -g -pg")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

# Glob all the source files in their respective directories
file(GLOB src_ninteg "src/*.f90")

# Create a shared library for ninteg
add_library(ninteg SHARED ${src_ninteg})
set_target_properties(ninteg PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_ninteg)

# Ensure the library is installed inside the python package folder
install(TARGETS ninteg DESTINATION ninteg)