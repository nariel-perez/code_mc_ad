!C----------------------------------------------------------------------------
!C-----------------------------------------------------------------------------
!******************************************************************************
!** FICHE F.14.  ALGORITHM TO HANDLE INDICES IN CONSTANT MU VT MONTE CARLO   **
!** This FORTRAN code is intended to illustrate points made in the text.     **
!** To our knowledge it works correctly.  However it is the responsibility of **
!** the user to test it, if it is to be used in a research application.      **
!******************************************************************************

!C    *******************************************************************
!C    ** INDEX-HANDLING IN GRAND CANONICAL MONTE CARLO SIMULATION.     **
!C    **                                                               **
!C    ** ROUTINES SUPPLIED:                                            **
!C    **                                                               **
!C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
!C    **    ADDS AN ATOM TO THE ARRAY LOCATE.                          **
!C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
!C    **    REMOVES AN ATOM FROM THE ARRAY LOCATE.                     **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                   NUMBER OF ATOMS BEFORE TRIAL      **
!C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS           **
!C    ** INTEGER IPULL               INDEX OF ATOM FOR REMOVAL         **
!C    ** INTEGER NLOC                POSITION OF N IN LOCATE           **
!C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING TRIAL      **
!C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
!C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF AN ATOM  **
!C    ** REAL    RX(NMAX), ETC.      POSITIONS OF ATOMS                **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** ROUTINE ADD IS CALLED AFTER A SUCCESSFUL TRIAL ADDITION.      **
!C    ** ROUTINE REMOVE IS CALLED AFTER A SUCCESSFUL TRIAL REMOVAL.    **
!C    ** THE ARRAY LOCATE IS UPDATED IN EACH CASE.                     **
!C    *******************************************************************
!C-------------------------------------------------------------------------
!C-------------------------------------------------------------------------

!C    *******************************************************************
!C    ** SUBROUTINE TO ADD AN ATOM TO THE ARRAY LOCATE.                **
!C    **                                                               **
!C    ** THERE ARE N ATOMS IN THE SIMULATION BEFORE THE NEW ADDITION   **
!C    *******************************************************************
!C    ******************************************************************


SUBROUTINE ADD(MOLKIND)
    USE AdsorbateInput, ONLY: NATOM
    USE SimulationData, ONLY: RX1, RY1, RZ1, LOCATE, RX, RY, RZ, ANX, ANGY, ANZ, EXNEW, EYNEW, EZNEW, N
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: MOLKIND

    ! Variables locales
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    INTEGER :: IPULL, I
    INTEGER :: INEW, NTRIAL
    INTEGER :: POSXIN, POSZIN, POSYIN

    ! Calcular posiciones iniciales... esto lo puedo comentar..
    POSXIN = INT(RX1(1) * 10000)
    POSYIN = INT(RY1(1) * 1000)
    POSZIN = INT(RZ1(1) * 100)

    ! Incrementar el número de moléculas de tipo MOLKIND
    NTRIAL = N(MOLKIND) + 1
    INEW = LOCATE(NTRIAL, MOLKIND)

    ! Si INEW es 0, asignar un nuevo número al átomo
    IF (INEW == 0) THEN
        LOCATE(NTRIAL, MOLKIND) = NTRIAL
        INEW = NTRIAL
    END IF

    ! Añadir las nuevas posiciones al array
    DO I = 1, NATOM(MOLKIND)
        RX(INEW, I, MOLKIND) = RX1(I)
        RY(INEW, I, MOLKIND) = RY1(I)
        RZ(INEW, I, MOLKIND) = RZ1(I)
    END DO

    ! Si la molécula tiene más de un átomo, guardar las orientaciones
    IF (NATOM(MOLKIND) > 1) THEN
        ANX(INEW, MOLKIND) = EXNEW
        ANGY(INEW, MOLKIND) = EYNEW
        ANZ(INEW, MOLKIND) = EZNEW
    END IF

    RETURN
END SUBROUTINE ADD
