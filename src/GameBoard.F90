MODULE game_board_mod

  USE kind_mod

  IMPLICIT NONE

  TYPE game_board_type
    PRIVATE
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: state_board
    INTEGER(ikind), DIMENSION(:,:), ALLOCATABLE :: neighbor_counts
  CONTAINS
    PROCEDURE, PASS(this_board) :: initialize
    PROCEDURE, PASS(this_board) :: count_cells
    PROCEDURE, PASS(this_board) :: accumulate_neighbors
    PROCEDURE, PASS(this_board) :: update_state
    PROCEDURE, PASS(this_board) :: print_state
  END TYPE

CONTAINS

  SUBROUTINE initialize(this_board)

    CLASS(game_board_type), INTENT(OUT) :: this_board

    INTEGER(ikind_large), PARAMETER :: sizex = 10000, sizey = 10000
    REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: cell_activity

    ALLOCATE(this_board%state_board(sizex,sizey),                               &
      this_board%neighbor_counts(sizex,sizey), cell_activity(sizex,sizey))
    CALL RANDOM_NUMBER(cell_activity)
    this_board%state_board = NINT(cell_activity)


  END SUBROUTINE initialize

  INTEGER(ikind_large) FUNCTION count_cells(this_board) RESULT(count)

    CLASS(game_board_type), INTENT(IN) :: this_board

    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%neighbor_counts, 1)
    ny = SIZE(this_board%neighbor_counts, 2)
    count = 0
!$omp parallel do reduction(+:count) shared(this_board)
    DO iy = 1, ny
  !$omp simd
      DO ix = 1, nx
        count = count + this_board%state_board(ix,iy)
      ENDDO
  !$omp end simd
    ENDDO
!$omp end parallel do

  END FUNCTION count_cells

  SUBROUTINE accumulate_neighbors(this_board)

    CLASS(game_board_type), INTENT(INOUT) :: this_board

    INTEGER(ikind_large) :: nx, ny, ix, iy

    nx = SIZE(this_board%neighbor_counts, 1)
    ny = SIZE(this_board%neighbor_counts, 2)
!$omp parallel shared(this_board)
  !$omp single
    this_board%neighbor_counts(1,1) = this_board%state_board(nx,1)              &
      + this_board%state_board(2,1) + this_board%state_board(1,ny)              &
      + this_board%state_board(1,2)
    !$omp simd
    DO ix = 2, nx-1
      this_board%neighbor_counts(ix, 1) = this_board%state_board(ix-1,1)        &
      + this_board%state_board(ix+1, 1) + this_board%state_board(ix,ny)         &
      + this_board%state_board(ix,2)
    ENDDO
    !$omp end simd
    this_board%neighbor_counts(nx,1) = this_board%state_board(nx-1,1)           &
      + this_board%state_board(1,1) + this_board%state_board(nx,ny)             &
      + this_board%state_board(nx,2)
  !$omp end single nowait
  !$omp do
    DO iy = 2, ny-1
      this_board%neighbor_counts(1,iy) = this_board%state_board(nx,iy)          &
        + this_board%state_board(2,iy) + this_board%state_board(ix,iy-1)        &
        + this_board%state_board(ix,iy+1)
    !$omp simd
      DO ix = 2, nx-1
        this_board%neighbor_counts(ix,iy) = this_board%state_board(ix-1,iy)     &
          + this_board%state_board(ix+1,iy) + this_board%state_board(ix,iy-1)   &
          + this_board%state_board(ix,iy+1)
      ENDDO
    !$omp end simd
      this_board%neighbor_counts(nx,iy) = this_board%state_board(nx-1,iy)       &
        + this_board%state_board(1,iy) + this_board%state_board(ix,iy-1)        &
        + this_board%state_board(ix,iy+1)
    ENDDO
  !$omp end do nowait
  !$omp single
    this_board%neighbor_counts(1,ny) = this_board%state_board(nx,ny)            &
      + this_board%state_board(2,ny) + this_board%state_board(1,ny-1)           &
      + this_board%state_board(1,1)
    !$omp simd
    DO ix = 2, nx-1
      this_board%neighbor_counts(ix, ny) = this_board%state_board(ix-1,ny)      &
        + this_board%state_board(ix+1, ny) + this_board%state_board(ix,ny-1)    &
        + this_board%state_board(ix,1)
    ENDDO
    !$omp end simd
    this_board%neighbor_counts(nx,ny) = this_board%state_board(nx-1,ny)         &
      + this_board%state_board(1,ny) + this_board%state_board(nx,ny-1)          &
      + this_board%state_board(nx,1)
  !$omp end single nowait
!$omp end parallel

  END SUBROUTINE accumulate_neighbors

  SUBROUTINE update_state(this_board)

    CLASS(game_board_type), INTENT(INOUT) :: this_board

!$omp parallel workshare shared(this_board)
    this_board%state_board = new_state(this_board%neighbor_counts,              &
      this_board%state_board)
!$omp end parallel workshare

  END SUBROUTINE update_state

  ELEMENTAL INTEGER(ikind) FUNCTION new_state(neighbor_count, current_state)

    INTEGER(ikind), INTENT(IN) :: neighbor_count
    INTEGER(ikind), INTENT(IN) :: current_state

    IF(current_state /= 0) THEN
      IF((neighbor_count == 2) .or. (neighbor_count == 3)) THEN
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

  SUBROUTINE print_state(this_board)

    CLASS(game_board_type), INTENT(IN) :: this_board

    INTEGER(ikind_large) :: row

    DO row = 1, SIZE(this_board%state_board, 2)
      WRITE(*,*) this_board%state_board(:,row)
    ENDDO

  END SUBROUTINE print_state

END MODULE