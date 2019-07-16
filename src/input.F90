MODULE input_mod

  USE kind_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_input, cell_pattern, sizex, sizey, nsteps

  CHARACTER(len=32) :: input_file_name = 'game_grid.nml'
  CHARACTER(len=32) :: cell_pattern = 'block', boundary = 'periodic'
  INTEGER(ikind_large) :: sizex = 3, sizey = 3, nsteps=3
  NAMELIST /starting_grid/ sizex, sizey, cell_pattern, nsteps, boundary

CONTAINS

  SUBROUTINE read_input()

    INTEGER(ikind) :: input_unit = 1
    LOGICAL :: file_exists

    IF(COMMAND_ARGUMENT_COUNT() > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, input_file_name)
    ENDIF
    INQUIRE(FILE=input_file_name, EXIST=file_exists)
    IF(file_exists) THEN
      OPEN(FILE=input_file_name, STATUS='OLD', NEWUNIT=input_unit)
      READ(UNIT=input_unit, NML=starting_grid)
    ENDIF

  END SUBROUTINE read_input

END MODULE input_mod