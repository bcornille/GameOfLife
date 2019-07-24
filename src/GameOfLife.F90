PROGRAM game_of_life

  USE kind_mod
  USE input_mod, ONLY: nsteps
  USE game_board_mod
#ifdef USE_MPI
  USE mpi_f08
#endif

  IMPLICIT NONE

  TYPE(game_board_type) :: game_board
  INTEGER(ikind_large) :: population, step

#ifdef USE_MPI
  DOUBLE PRECISION :: tstart, total_time
  INTEGER :: ierror, mpi_thread_provided, rank
  CALL MPI_Init_thread(MPI_THREAD_MULTIPLE, mpi_thread_provided, ierror)
  IF(mpi_thread_provided < MPI_THREAD_FUNNELED) THEN
    CALL MPI_Finalize(ierror)
    STOP "Not a threaded MPI implementation."
  ENDIF
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
#else
  INTEGER(ikind_large) :: count_start, count_end, count_rate
  REAL(rkind) :: total_time
#endif

  !Start by initializing the game board.
  CALL game_board%initialize()
  population = game_board%count_cells()
  ! CALL game_board%print_state()
#ifdef USE_MPI
  IF(rank == 0) THEN
#endif
  WRITE(*,'(A12,I8)') "Population: ", population
#ifdef USE_MPI
  ENDIF
  tstart = MPI_Wtime()
#else
  CALL SYSTEM_CLOCK(count_start, count_rate)
#endif

  DO step = 1, nsteps
    CALL game_board%update_halo()
    ! CALL game_board%print_state()
    CALL game_board%accumulate_neighbors()
    CALL game_board%update_state()
    population = game_board%count_cells()
    ! CALL game_board%print_state()
#ifdef USE_MPI
    IF(rank == 0) THEN
#endif
    WRITE(*,'(A12,I8)') "Population: ", population
#ifdef USE_MPI
    ENDIF
#endif
  ENDDO
#ifdef USE_MPI
  total_time = MPI_Wtime() - tstart
#else
  CALL SYSTEM_CLOCK(count_end)
  total_time = (count_end - count_start) / REAL(count_rate, rkind)
#endif
  CALL game_board%update_halo()
  ! CALL game_board%print_state()
#ifdef USE_MPI
  IF(rank == 0) THEN
#endif
  WRITE(*,'(A9,F8.3,A8)') "Runtime: ", total_time, " seconds"
#ifdef USE_MPI
  ENDIF
#endif

#ifdef USE_MPI
  CALL MPI_Finalize(ierror)
#endif

END PROGRAM game_of_life
