MODULE kind_mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: size_ikind = 1,                                         &
                        ikind = SELECTED_INT_KIND(size_ikind),                  &
                        ikind_large = SELECTED_INT_KIND(8),                     &
                        rkind = SELECTED_REAL_KIND(15,1)

END MODULE kind_mod
