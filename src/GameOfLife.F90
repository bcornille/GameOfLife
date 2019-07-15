PROGRAM game_of_life

  USE game_board_mod

  IMPLICIT NONE

  TYPE(game_board_type) :: game_board
  INTEGER(ikind_large) :: population, nsteps = 50, step

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

END PROGRAM game_of_life