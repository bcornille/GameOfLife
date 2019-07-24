# GameOfLife
Conway's "Game of Life" Implemented in Modern Fortran with MPI and OpenMP

## Features Demonstrated

### OO-Fortran

The `game_board_type` class contains all the necessary data and all operations
on these data are type-bound procedures. The overall program should be self-
documenting as the basic steps of the algorithm are appropriately named
procedures (although the code is somewhat obfuscated by MPI preprocessor
guards).

### MPI

This code uses the `mpi_f08` interface which has a more mondern Fortran style.
It also should allow for the use of passing array sections to MPI routines
without the cost of unnecessary copies. Support for this style appears to be
very recent with some compilers. Some of the MPI features of note that are used
include:

  - Thread-aware MPI initialization
  - Fortran 90 kind equivalent `MPI_Datatype` creation
  - Creation of a Cartesian MPI communicator to allow more optimal rank ordering
    for the network topology (performance impact unmeasured)
  - Non-blocking send/receive calls
  - Collective reduction

### OpenMP

Some intermediate features of OpenMP are used.

  - `collapse` of `parallel do` loops for nested `DO` blocks
  - conditional entry into OpenMP regions
  - `parallel sections` for work-sharing of a set number of independent
    operations; used for MPI send/receive calls when the MPI implementation
    supports calls from multiple threads (could be performance gain or loss)
  - `reduction` clause in `parallel do` statement
