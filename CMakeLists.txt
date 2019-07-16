cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(GameOfLife VERSION 0.1 LANGUAGES Fortran)

set(ENABLE_OpenMP FALSE CACHE BOOL "Whether to enable OpenMP: Default is false")
set(ENABLE_MPI FALSE CACHE BOOL "Whether to enable OpenMP: Default is false")

if(ENABLE_OpenMP)
	find_package(OpenMP)
endif()
if(ENABLE_MPI)
	find_package(MPI)
endif()

add_subdirectory(src)