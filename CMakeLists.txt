cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(GameOfLife VERSION 0.1 LANGUAGES Fortran)

set(ENABLE_OpenMP FALSE CACHE BOOL "Whether to enable OpenMP: Default is false")
set(ENABLE_MPI FALSE CACHE BOOL "Whether to enable OpenMP: Default is false")

if(ENABLE_OpenMP)
	find_package(OpenMP)
endif()
if(ENABLE_MPI)
	find_package(MPI REQUIRED)
	if(MPI_Fortran_FOUND)
		add_definitions(-DUSE_MPI)
	endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")

add_subdirectory(src)