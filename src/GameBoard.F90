MODULE game_board_mod

  USE kind_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: game_board_type

  TYPE game_board_type
    PRIVATE
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: state_board
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: neighbor_counts
  CONTAINS
    PROCEDURE, PASS(this_board) :: initialize
    PROCEDURE, PASS(this_board) :: count_cells
    PROCEDURE, PASS(this_board) :: accumulate_neighbors
    PROCEDURE, PASS(this_board) :: update_state
    PROCEDURE, PASS(this_board), PRIVATE :: update_halo
    PROCEDURE, PASS(this_board) :: print_state
  END TYPE

CONTAINS

  SUBROUTINE initialize(this_board)

    USE input_mod, ONLY: read_input, cell_pattern, sizex, sizey

    CLASS(game_board_type), INTENT(OUT) :: this_board

    INTEGER(ikind_large) :: ix, iy
    REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: cell_activity

    CALL read_input()
    ALLOCATE(this_board%state_board(sizex+2,sizey+2),                           &
          this_board%neighbor_counts(sizex,sizey))
    this_board%state_board = 0
    SELECT CASE(cell_pattern)
      CASE('block')
        IF((sizex < 3) .OR. (sizey < 3)) THEN
          STOP "Grid size too small for cell pattern."
        ENDIF
        this_board%state_board(1:3,1:3) = RESHAPE(                              &
          [[0, 0, 0],                                                           &
           [0, 1, 1],                                                           &
           [0, 1, 1]], [3,3])
      CASE('blinker')
        IF((sizex < 4) .OR. (sizey < 4)) THEN
          STOP "Grid size too small for cell pattern."
        ENDIF
        this_board%state_board(1:4,1:4) = RESHAPE(                              &
          [[0, 0, 0, 0],                                                        &
           [0, 0, 1, 0],                                                        &
           [0, 0, 1, 0],                                                        &
           [0, 0, 1, 0]], [4,4])
      CASE('glider')
        IF((sizex < 5) .OR. (sizey < 5)) THEN
          STOP "Grid size too small for cell pattern."
        ENDIF
        this_board%state_board(1:4,1:4) = RESHAPE(                              &
          [[0, 0, 0, 0],                                                        &
           [0, 0, 1, 0],                                                        &
           [0, 0, 0, 1],                                                        &
           [0, 1, 1, 1]], [4,4])
      CASE('R-pentomino')
        IF((sizex < 5) .OR. (sizey < 5)) THEN
          STOP "Grid size too small for cell pattern."
        ENDIF
        this_board%state_board(1:4,1:4) = RESHAPE(                              &
          [[0, 0, 0, 0],                                                        &
           [0, 0, 1, 1],                                                        &
           [0, 1, 1, 0],                                                        &
           [0, 0, 1, 0]], [4,4])

      CASE('die-hard')
        IF((sizex < 10) .OR. (sizey < 5)) THEN
          STOP "Grid size too small for cell pattern."
        ENDIF
        this_board%state_board(1:9,1:4) = RESHAPE(                              &
          [[0, 0, 0, 0, 0, 0, 0, 0, 0],                                                        &
           [0, 0, 0, 0, 0, 0, 0, 1, 0],                                                        &
           [0, 1, 1, 0, 0, 0, 0, 0, 0],                                                        &
           [0, 0, 1, 0, 0, 0, 1, 1, 1]], [9,4])
      CASE('random')
        ALLOCATE(cell_activity(sizex,sizey))
!$omp parallel do collapse(2) shared(this_board, cell_activity)
        DO iy = 2, sizey-1
          DO ix = 2, sizex-1
            CALL RANDOM_NUMBER(cell_activity(ix-1,iy-1))
            this_board%state_board(ix,iy) = NINT(cell_activity(ix-1,iy-1))
          ENDDO
        ENDDO
!$omp end parallel do
      CASE DEFAULT
        STOP "Cell pattern not recognized."
    END SELECT

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

  END FUNCTION count_cells

  SUBROUTINE accumulate_neighbors(this_board)

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)
    ASSOCIATE (state => this_board%state_board)
!$omp parallel do collapse(2) shared(this_board)
      DO iy = 2, ny-1
        DO ix = 2, nx-1
          this_board%neighbor_counts(ix-1,iy-1) =                                   &
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
    CALL this_board%update_halo()

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

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER(ikind_large) :: nx, ny

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)

!$omp parallel shared(this_board)
  !$omp single
    this_board%state_board(1,1) = this_board%state_board(nx-1,ny-1)
  !$omp end single nowait
  !$omp single
    this_board%state_board(2:nx-1,1) = this_board%state_board(2:nx-1,ny-1)
  !$omp end single nowait
  !$omp single
    this_board%state_board(nx,1) = this_board%state_board(2,ny-1)
  !$omp end single nowait
  !$omp single
    this_board%state_board(1,2:ny-1) = this_board%state_board(nx-1,2:ny-1)
  !$omp end single nowait
  !$omp single
    this_board%state_board(nx,2:ny-1) = this_board%state_board(2,2:ny-1)
  !$omp end single nowait
  !$omp single
    this_board%state_board(1,ny) = this_board%state_board(nx-1,2)
  !$omp end single nowait
  !$omp single
    this_board%state_board(2:nx-1,ny) = this_board%state_board(2:nx-1,2)
  !$omp end single nowait
  !$omp single
    this_board%state_board(nx,ny) = this_board%state_board(2,2)
  !$omp end single nowait
!$omp end parallel

  END SUBROUTINE update_halo

  SUBROUTINE print_state(this_board)

    CLASS(game_board_type), INTENT(IN) :: this_board

    INTEGER(ikind_large) :: nx, ny, row

    nx = SIZE(this_board%state_board, 1)
    ny = SIZE(this_board%state_board, 2)

    DO row = 2, ny-1
      WRITE(*,*) this_board%state_board(2:nx-1,row)
    ENDDO

  END SUBROUTINE print_state

END MODULE