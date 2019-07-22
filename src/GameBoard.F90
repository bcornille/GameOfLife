MODULE game_board_mod

  USE kind_mod
#ifdef USE_MPI
  USE mpi_f08
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: game_board_type

#ifdef USE_MPI
  TYPE mpi_grid_type
    INTEGER :: my_rank, my_coord(2), num_ranks, provided
    TYPE(MPI_Request) :: send_requests(8), recv_requests(8)
    INTEGER :: northwest, north, northeast, west,                               &
      east, southwest, south, southeast
    TYPE(MPI_Comm) :: comm
    TYPE(MPI_Datatype) :: corner_type
  CONTAINS
    PROCEDURE, PASS(this_grid) :: setup_comm
  END TYPE mpi_grid_type
#endif

  TYPE game_board_type
    PRIVATE
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: state_board
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: neighbor_counts
#ifdef USE_MPI
    TYPE(mpi_grid_type), PUBLIC :: mpi_grid
#endif
  CONTAINS
    PROCEDURE, PASS(this_board) :: initialize
    PROCEDURE, PASS(this_board) :: count_cells
    PROCEDURE, PASS(this_board) :: accumulate_neighbors
    PROCEDURE, PASS(this_board) :: update_state
    PROCEDURE, PASS(this_board) :: update_halo
    PROCEDURE, PASS(this_board) :: print_state
  END TYPE

CONTAINS

  SUBROUTINE initialize(this_board)

    USE input_mod, ONLY: read_input, cell_pattern, sizex, sizey,                &
      offsetx, offsety, init

    CLASS(game_board_type), INTENT(OUT) :: this_board

    INTEGER :: ierror
    INTEGER(ikind_large) :: ix, iy
    REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: cell_activity

    CALL read_input()
#ifdef USE_MPI
    CALL this_board%mpi_grid%setup_comm()
#endif
    ALLOCATE(this_board%state_board(sizex+2,sizey+2),                           &
          this_board%neighbor_counts(sizex,sizey))
    this_board%state_board = 0
#ifdef USE_MPI
    IF((this_board%mpi_grid%my_rank == 0).OR.(TRIM(init) == 'all')) THEN
#endif
    SELECT CASE(cell_pattern)
      CASE('block')
        IF((sizex < 3 + offsetx) .OR. (sizey < 3 + offsety)) THEN
#ifdef USE_MPI
          CALL MPI_Finalize()
#endif
          STOP "Grid size too small for cell pattern "//TRIM(cell_pattern)//"."
        ENDIF
        this_board%state_board(1+offsetx:3+offsetx,1+offsety:3+offsety) =       &
        RESHAPE([[0, 0, 0],                                                     &
                 [0, 1, 1],                                                     &
                 [0, 1, 1]], [3,3])
      CASE('blinker')
        IF((sizex < 4 + offsetx) .OR. (sizey < 4 + offsety)) THEN
#ifdef USE_MPI
          CALL MPI_Finalize()
#endif
          STOP "Grid size too small for cell pattern "//TRIM(cell_pattern)//"."
        ENDIF
        this_board%state_board(1+offsetx:4+offsetx,1+offsety:4+offsety) =       &
        RESHAPE([[0, 0, 0, 0],                                                  &
                 [0, 0, 1, 0],                                                  &
                 [0, 0, 1, 0],                                                  &
                 [0, 0, 1, 0]], [4,4])
      CASE('glider')
        IF((sizex < 5 + offsetx) .OR. (sizey < 5 + offsety)) THEN
#ifdef USE_MPI
          CALL MPI_Finalize()
#endif
          STOP "Grid size too small for cell pattern "//TRIM(cell_pattern)//"."
        ENDIF
        this_board%state_board(1+offsetx:4+offsetx,1+offsety:4+offsety) =       &
        RESHAPE([[0, 0, 0, 0],                                                  &
                 [0, 0, 1, 0],                                                  &
                 [0, 0, 0, 1],                                                  &
                 [0, 1, 1, 1]], [4,4])
      CASE('R-pentomino')
        IF((sizex < 5 + offsetx) .OR. (sizey < 5 + offsety)) THEN
#ifdef USE_MPI
          CALL MPI_Finalize()
#endif
          STOP "Grid size too small for cell pattern "//TRIM(cell_pattern)//"."
        ENDIF
        this_board%state_board(1+offsetx:4+offsetx,1+offsety:4+offsety) =       &
        RESHAPE([[0, 0, 0, 0],                                                  &
                 [0, 0, 1, 1],                                                  &
                 [0, 1, 1, 0],                                                  &
                 [0, 0, 1, 0]], [4,4])

      CASE('die-hard')
        IF((sizex < 10 + offsetx) .OR. (sizey < 5 + offsety)) THEN
#ifdef USE_MPI
          CALL MPI_Finalize()
#endif
          STOP "Grid size too small for cell pattern "//TRIM(cell_pattern)//"."
        ENDIF
        this_board%state_board(1+offsetx:9+offsetx,1+offsety:4+offsety) =       &
        RESHAPE([[0, 0, 0, 0, 0, 0, 0, 0, 0],                                   &
                 [0, 0, 0, 0, 0, 0, 0, 1, 0],                                   &
                 [0, 1, 1, 0, 0, 0, 0, 0, 0],                                   &
                 [0, 0, 1, 0, 0, 0, 1, 1, 1]], [9,4])
      CASE('random')
        ALLOCATE(cell_activity(sizex,sizey))
!$omp parallel do collapse(2) shared(this_board, cell_activity)
        DO iy = 1, sizey
          DO ix = 1, sizex
            CALL RANDOM_NUMBER(cell_activity(ix,iy))
            this_board%state_board(ix+1,iy+1) = NINT(cell_activity(ix,iy))
          ENDDO
        ENDDO
!$omp end parallel do
      CASE DEFAULT
#ifdef USE_MPI
        CALL MPI_Finalize()
#endif
        STOP "Cell pattern "//TRIM(cell_pattern)//" not recognized."
    END SELECT
#ifdef USE_MPI
    ENDIF
    ASSOCIATE(state => this_board%state_board, grid => this_board%mpi_grid)
!$omp parallel sections if(grid%provided == MPI_THREAD_MULTIPLE) shared(this_board)
 !$omp section
      CALL MPI_ISend(state(2,2), 1, grid%corner_type, grid%northwest, 0,        &
        grid%comm, grid%send_requests(1), ierror)
 !$omp section
      CALL MPI_ISend(state(2:sizex+1,2), sizex, grid%corner_type, grid%north,  &
        1, grid%comm, grid%send_requests(2), ierror)
 !$omp section
      CALL MPI_ISend(state(sizex+1,2), 1, grid%corner_type, grid%northeast, 2,  &
        grid%comm, grid%send_requests(3), ierror)
 !$omp section
      CALL MPI_ISend(state(2,2:sizey+1), sizey, grid%corner_type, grid%west, 3,  &
        grid%comm, grid%send_requests(4), ierror)
 !$omp section
      CALL MPI_ISend(state(sizex+1,2:sizey+1), sizey, grid%corner_type,          &
        grid%east, 4, grid%comm, grid%send_requests(5), ierror)
 !$omp section
      CALL MPI_ISend(state(2,sizey+1), 1, grid%corner_type, grid%southwest, 5,  &
        grid%comm, grid%send_requests(6), ierror)
 !$omp section
      CALL MPI_ISend(state(2:sizex+1,sizey+1), sizex, grid%corner_type,        &
        grid%south, 6, grid%comm, grid%send_requests(7), ierror)
 !$omp section
      CALL MPI_ISend(state(sizex+1,sizey+1), 1, grid%corner_type,               &
        grid%southeast, 7, grid%comm, grid%send_requests(8), ierror)
!$omp end parallel sections
    END ASSOCIATE
#endif

  END SUBROUTINE initialize

  INTEGER(ikind_large) FUNCTION count_cells(this_board) RESULT(count)

    CLASS(game_board_type), INTENT(IN) :: this_board

    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)
    count = 0
!$omp parallel do collapse(2) reduction(+:count) shared(this_board)
    DO iy = 2, ny-1
      DO ix = 2, nx-1
        count = count + this_board%state_board(ix,iy)
      ENDDO
    ENDDO
!$omp end parallel do
#ifdef USE_MPI
    CALL MPI_Allreduce(MPI_IN_PLACE, count, 1, MPI_INTEGER, MPI_SUM,            &
      this_board%mpi_grid%comm)
#endif

  END FUNCTION count_cells

  SUBROUTINE accumulate_neighbors(this_board)

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)
    ASSOCIATE(state => this_board%state_board)
!$omp parallel do collapse(2) shared(this_board)
      DO iy = 2, ny-1
        DO ix = 2, nx-1
          this_board%neighbor_counts(ix-1,iy-1) =                               &
            state(ix-1,iy-1) + state(ix,iy-1) + state(ix+1,iy-1)                &
          + state(ix-1,iy)                    + state(ix+1,iy)                  &
          + state(ix-1,iy+1) + state(ix,iy+1) + state(ix+1,iy+1)
        ENDDO
      ENDDO
!$omp end parallel do
    END ASSOCIATE

  END SUBROUTINE accumulate_neighbors

  SUBROUTINE update_state(this_board)

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER :: ierror
    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)

!$omp parallel do collapse(2) shared(this_board)
    DO iy = 2, ny-1
      DO ix = 2, nx-1
        this_board%state_board(ix,iy) = new_state(                              &
          this_board%neighbor_counts(ix-1,iy-1), this_board%state_board(ix,iy))
      ENDDO
    ENDDO
!$omp end parallel do
#ifdef USE_MPI
    ASSOCIATE(state => this_board%state_board, grid => this_board%mpi_grid)
!$omp parallel sections if(grid%provided == MPI_THREAD_MULTIPLE) shared(this_board)
 !$omp section
      CALL MPI_ISend(state(2,2), 1, grid%corner_type, grid%northwest, 0,        &
        grid%comm, grid%send_requests(1), ierror)
 !$omp section
      CALL MPI_ISend(state(2:nx-1,2), nx-2, grid%corner_type, grid%north, 1,  &
        grid%comm, grid%send_requests(2), ierror)
 !$omp section
      CALL MPI_ISend(state(nx-1,2), 1, grid%corner_type, grid%northeast, 2,     &
        grid%comm, grid%send_requests(3), ierror)
 !$omp section
      CALL MPI_ISend(state(2,2:ny-1), ny-2, grid%corner_type, grid%west, 3,     &
        grid%comm, grid%send_requests(4), ierror)
 !$omp section
      CALL MPI_ISend(state(nx-1,2:ny-1), ny-2, grid%corner_type, grid%east, 4,  &
        grid%comm, grid%send_requests(5), ierror)
 !$omp section
      CALL MPI_ISend(state(2,ny-1), 1, grid%corner_type, grid%southwest, 5,     &
        grid%comm, grid%send_requests(6), ierror)
 !$omp section
      CALL MPI_ISend(state(2:nx-1,ny-1), nx-2, grid%corner_type, grid%south,  &
        6, grid%comm, grid%send_requests(7), ierror)
 !$omp section
      CALL MPI_ISend(state(nx-1,ny-1), 1, grid%corner_type, grid%southeast, 7,  &
        grid%comm, grid%send_requests(8), ierror)
!$omp end parallel sections
    END ASSOCIATE
#endif

  END SUBROUTINE update_state

  ELEMENTAL INTEGER(ikind) FUNCTION new_state(neighbor_count, current_state)

    INTEGER(ikind), INTENT(IN) :: neighbor_count
    INTEGER(ikind), INTENT(IN) :: current_state

    IF(current_state /= 0) THEN
      IF((neighbor_count == 2) .OR. (neighbor_count == 3)) THEN
        new_state = 1
      ELSE
        new_state = 0
      ENDIF
    ELSE
      IF(neighbor_count == 3) THEN
        new_state = 1
      ELSE
        new_state = 0
      ENDIF
    ENDIF

  END FUNCTION new_state

  SUBROUTINE update_halo(this_board)

    USE input_mod, ONLY: boundary

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER(ikind_large) :: nx, ny, ierror

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)

#ifndef USE_MPI
    SELECT CASE(boundary)
      CASE ('periodic')
!!$omp parallel sections shared(this_board)
!  !$omp section
        this_board%state_board(1,1) = this_board%state_board(nx-1,ny-1)
!  !$omp section
        this_board%state_board(2:nx-1,1) = this_board%state_board(2:nx-1,ny-1)
!  !$omp section
        this_board%state_board(nx,1) = this_board%state_board(2,ny-1)
!  !$omp section
        this_board%state_board(1,2:ny-1) = this_board%state_board(nx-1,2:ny-1)
!  !$omp section
        this_board%state_board(nx,2:ny-1) = this_board%state_board(2,2:ny-1)
!  !$omp section
        this_board%state_board(1,ny) = this_board%state_board(nx-1,2)
!  !$omp section
        this_board%state_board(2:nx-1,ny) = this_board%state_board(2:nx-1,2)
!  !$omp section
        this_board%state_board(nx,ny) = this_board%state_board(2,2)
!!$omp end parallel sections
      CASE('dead')
      CASE DEFAULT
        STOP "Boundary condition "//TRIM(boundary)//" not implemented."
    END SELECT
#else
    ASSOCIATE(state => this_board%state_board, grid => this_board%mpi_grid)
!$omp parallel sections if(grid%provided == MPI_THREAD_MULTIPLE)
 !$omp section
      CALL MPI_IRecv(state(1,1), 1, grid%corner_type, grid%northwest, 7,         &
        grid%comm, grid%recv_requests(1), ierror)
 !$omp section
      CALL MPI_IRecv(state(2:nx-1,1), nx+2, grid%corner_type, grid%north, 6,   &
        grid%comm, grid%recv_requests(2), ierror)
 !$omp section
      CALL MPI_IRecv(state(nx,1), 1, grid%corner_type, grid%northeast, 5,        &
        grid%comm, grid%recv_requests(3), ierror)
 !$omp section
      CALL MPI_IRecv(state(1,2:ny-1), ny+2, grid%corner_type, grid%west, 4,      &
        grid%comm, grid%recv_requests(4), ierror)
 !$omp section
      CALL MPI_IRecv(state(nx,2:ny-1), ny+2, grid%corner_type, grid%east, 3,     &
        grid%comm, grid%recv_requests(5), ierror)
 !$omp section
      CALL MPI_IRecv(state(1,ny), 1, grid%corner_type, grid%southwest, 2,        &
        grid%comm, grid%recv_requests(6), ierror)
 !$omp section
      CALL MPI_IRecv(state(2:nx-1,ny), nx+2, grid%corner_type, grid%south, 1,  &
        grid%comm, grid%recv_requests(7), ierror)
 !$omp section
      CALL MPI_IRecv(state(nx,ny), 1, grid%corner_type, grid%southeast, 0,        &
        grid%comm, grid%recv_requests(8), ierror)
 !$omp section
      CALL MPI_Waitall(8, grid%send_requests, MPI_STATUSES_IGNORE)
!$omp end parallel sections
      CALL MPI_Waitall(8, grid%recv_requests, MPI_STATUSES_IGNORE)
    END ASSOCIATE
#endif

  END SUBROUTINE update_halo

  SUBROUTINE print_state(this_board)

    CLASS(game_board_type), INTENT(IN) :: this_board

    INTEGER(ikind_large) :: nx, ny, row
    INTEGER :: rank, ierror

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)

#ifdef USE_MPI
    DO rank = 0, this_board%mpi_grid%num_ranks - 1
      IF(this_board%mpi_grid%my_rank == rank) THEN
#endif
    DO row = 1, ny
      WRITE(*,*) this_board%state_board(:,row)
    ENDDO
    WRITE(*,*)
#ifdef USE_MPI
      ENDIF
      CALL MPI_Barrier(this_board%mpi_grid%comm, ierror)
    ENDDO
    CALL MPI_Barrier(this_board%mpi_grid%comm, ierror)
#endif

  END SUBROUTINE print_state

#ifdef USE_MPI
  SUBROUTINE setup_comm(this_grid)

    USE input_mod, ONLY: boundary, sizex, sizey

    CLASS(mpi_grid_type) :: this_grid

    LOGICAL :: periodic(2) = .FALSE., reorder = .TRUE.
    INTEGER :: ierror, dims(2) = 0
    INTEGER(ikind_large) :: remainder

    CALL MPI_Query_thread(this_grid%provided, ierror)
    CALL MPI_Type_create_f90_integer(1, this_grid%corner_type, ierror)
    CALL MPI_Type_commit(this_grid%corner_type, ierror)
    IF(TRIM(boundary) == 'periodic') THEN
      periodic = .TRUE.
    ENDIF
    CALL MPI_Comm_size(MPI_COMM_WORLD, this_grid%num_ranks, ierror)
    CALL MPI_Dims_create(this_grid%num_ranks, 2, dims, ierror)
    CALL MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, reorder,            &
      this_grid%comm, ierror)
    CALL MPI_Comm_rank(this_grid%comm, this_grid%my_rank, ierror)
    this_grid%my_coord(1) = MOD(this_grid%my_rank, dims(1))
    this_grid%my_coord(2) = this_grid%my_rank / dims(1)
    remainder = MOD(sizex, dims(1))
    sizex = sizex / dims(1)
    IF(this_grid%my_coord(1) < remainder) THEN
      sizex = sizex + 1
    ENDIF
    remainder = MOD(sizey, dims(2))
    sizey = sizey / dims(2)
    IF(this_grid%my_coord(2) < remainder) THEN
      sizey = sizey + 1
    ENDIF
    ! WRITE(*,*) this_grid%my_rank, sizex, sizey
    this_grid%northwest = MOD(this_grid%my_coord(1) + dims(1) - 1, dims(1))     &
      + MOD(this_grid%my_coord(2) + dims(2) - 1, dims(2))*dims(1)
    this_grid%north = this_grid%my_coord(1)                                     &
      + MOD(this_grid%my_coord(2) + dims(2) - 1, dims(2))*dims(1)
    this_grid%northeast = MOD(this_grid%my_coord(1) + dims(1) + 1, dims(1))     &
      + MOD(this_grid%my_coord(2) + dims(2) - 1, dims(2))*dims(1)
    this_grid%west = this_grid%my_coord(2)*dims(1)                              &
      + MOD(this_grid%my_coord(1) + dims(1) - 1, dims(1))
    this_grid%east = this_grid%my_coord(2)*dims(1)                              &
      + MOD(this_grid%my_coord(1) + dims(1) + 1, dims(1))
    this_grid%southwest = MOD(this_grid%my_coord(1) + dims(1) - 1, dims(1))     &
      + MOD(this_grid%my_coord(2) + dims(2) + 1, dims(2))*dims(1)
    this_grid%south = this_grid%my_coord(1)                                     &
      + MOD(this_grid%my_coord(2) + dims(2) + 1, dims(2))*dims(1)
    this_grid%southeast = MOD(this_grid%my_coord(1) + dims(1) + 1, dims(1))     &
      + MOD(this_grid%my_coord(2) + dims(2) + 1, dims(2))*dims(1)
    IF(TRIM(boundary) /= 'periodic') THEN
      IF(this_grid%my_coord(1) == 0) THEN
        this_grid%northwest = MPI_PROC_NULL
        this_grid%west = MPI_PROC_NULL
        this_grid%southwest = MPI_PROC_NULL
      ENDIF
      IF(this_grid%my_coord(1) == dims(1) - 1) THEN
        this_grid%northeast = MPI_PROC_NULL
        this_grid%east = MPI_PROC_NULL
        this_grid%southeast = MPI_PROC_NULL
      ENDIF
      IF(this_grid%my_coord(2) == 0) THEN
        this_grid%northwest = MPI_PROC_NULL
        this_grid%north = MPI_PROC_NULL
        this_grid%northeast = MPI_PROC_NULL
      ENDIF
      IF(this_grid%my_coord(2) == dims(2) - 1) THEN
        this_grid%southwest = MPI_PROC_NULL
        this_grid%south = MPI_PROC_NULL
        this_grid%southeast = MPI_PROC_NULL
      ENDIF
    ENDIF
    ! IF(this_grid%my_rank == 0) THEN
    !   WRITE(*,*) this_grid%northwest, this_grid%north, this_grid%northeast
    !   WRITE(*,*) this_grid%west, this_grid%east
    !   WRITE(*,*) this_grid%southwest, this_grid%south, this_grid%southeast
    ! ENDIF

  END SUBROUTINE setup_comm
#endif

END MODULE