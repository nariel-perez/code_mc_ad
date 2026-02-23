REAL FUNCTION RANF(DUMMY)
    IMPLICIT NONE
    REAL :: DUMMY

    !*******************************************************************
    ! ** RETURNS A UNIFORM RANDOM VARIATE IN THE RANGE 0 TO 1.         **
    ! **                                                            **
    ! **                 ***************                               **
    ! **                 **  WARNING  **                               **
    ! **                 ***************                               **
    ! **                                                               **
    ! ** GOOD RANDOM NUMBER GENERATORS ARE MACHINE SPECIFIC.           **
    ! ** PLEASE USE THE ONE RECOMMENDED FOR YOUR MACHINE.              **
    !*******************************************************************

    CALL RANDOM_NUMBER(RANF)

    RETURN
END FUNCTION RANF
