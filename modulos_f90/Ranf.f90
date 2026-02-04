REAL FUNCTION RANF(DUMMY)
    IMPLICIT NONE
    ! Parámetros para el generador de números aleatorios
    INTEGER :: L, C, M
    PARAMETER (L = 1029, C = 221591, M = 1048576)

    ! Variables locales
    INTEGER :: SEED
    REAL :: DUMMY
    SAVE :: SEED
    DATA SEED / 0 /

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

    ! Actualizar la semilla usando el método lineal congruencial
    SEED = MOD(SEED * L + C, M)

    ! Calcular el número aleatorio en el rango [0, 1)
    RANF = REAL(SEED) / M

    RETURN
END FUNCTION RANF
