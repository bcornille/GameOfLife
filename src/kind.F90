MODULE kind_mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(1),                            &
                        ikind_large = SELECTED_INT_KIND(8),                     &
                        rkind = SELECTED_REAL_KIND(15,1)

END MODULE kind_mod