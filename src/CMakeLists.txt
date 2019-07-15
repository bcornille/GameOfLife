set(GOL_SOURCES
	GameOfLife.F90
	GameBoard.F90
	kind.F90
)

add_executable(GameOfLife ${GOL_SOURCES})
if(OpenMP_FOUND)
	target_link_libraries(GameOfLife OpenMP::OpenMP_Fortran)
endif()