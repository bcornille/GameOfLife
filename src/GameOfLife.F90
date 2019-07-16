PROGRAM game_of_life

  USE kind_mod
  USE input_mod, ONLY: nsteps
  USE game_board_mod

  IMPLICIT NONE

  TYPE(game_board_type) :: game_board
  INTEGER(ikind_large) :: population, step

#if USE_MPI
  INTEGER :: mpi_error
  MPI_Init(mpi_error)
#endif

  !Start by initializing the game board.
  CALL game_board%initialize()
  population = game_board%count_cells()
  ! CALL game_board%print_state()
  WRITE(*,*) "Population: ", population


  DO step = 1, nsteps
    CALL game_board%accumulate_neighbors()
    CALL game_board%update_state()
    population = game_board%count_cells()
    ! CALL game_board%print_state()
    WRITE(*,*) "Population: ", population
  ENDDO

#if USE_MPI
  MPI_Finalize(mpi_error)
#endif

END PROGRAM game_of_life