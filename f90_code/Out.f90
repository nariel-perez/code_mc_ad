!C----------------------------------------------------------------------
!C-----------------------------------------------------------------------
!C    *******************************************************************
!C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** REAL    TEMP         TEMPERATURE                              **
!C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
!C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
!C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
!C    ** REAL    V            POTENTIAL ENERGY                         **
!C    ** REAL    W            VIRIAL                                   **
!C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
!C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
!C    *******************************************************************
!----------------------------------------------------


SUBROUTINE OUT(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, GHOST, JPASOS, CANONICALMOLECULES)

    USE SimulationData, ONLY: NMOLEC, N, LOCATE
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: JPASOS, CANONICALMOLECULES
    REAL :: TEMP, Z(10), SIGMA, EPS, RCUT
    REAL :: V, VA, VG, W
    LOGICAL :: GHOST

    ! Variables locales
    INTEGER :: MOLKIND
    INTEGER :: NCONFMIN, NCONFMAX
    INTEGER :: NTOTALB, I
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    REAL :: BETA, DELTV, DELTW, DELTDB, DELTVA
    REAL :: RANF, DUMMY
    INTEGER :: NTRIAL, NLOC, IPULL
    LOGICAL :: OVRLAP, CREATE
    REAL :: sorteo, b
    REAL :: Ntotalmolec
    REAL :: ZTOTAL

    ! Inicialización de variables
    CREATE = .FALSE.
    GHOST = .FALSE.
    BETA = 1.0 / TEMP

    ! Elegir tipo de molécula
    MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1
    NTRIAL = N(MOLKIND) - 1

    ! Seleccionar un elemento aleatorio de la parte activa de LOCATE
    b = RANF(DUMMY)
    NLOC = INT(REAL(NTRIAL) * b) + 1
    IPULL = LOCATE(NLOC, MOLKIND)

    ! Calcular el cambio de energía al eliminar el átomo IPULL
    CALL POTOUT(IPULL, MOLKIND, DELTV)
    CALL ADPOTOUT(IPULL, MOLKIND, DELTVA)

    DELTDB = BETA * (DELTV + DELTVA) - LOG(N(MOLKIND) / Z(MOLKIND))

    ! Verificar la aceptación del cambio
    IF (DELTDB < 75.0) THEN
        IF (DELTDB < 0.0) THEN
            GHOST = .TRUE.
            CALL REMOVE(NLOC, IPULL, MOLKIND)
            V = V + DELTV + DELTVA
            VG = VG + DELTV
            VA = VA + DELTVA
            N(MOLKIND) = NTRIAL
        ELSE
            sorteo = RANF(DUMMY)
            IF (EXP(-DELTDB) > sorteo) THEN
                GHOST = .TRUE.
                CALL REMOVE(NLOC, IPULL, MOLKIND)
                V = V + DELTV + DELTVA
                VG = VG + DELTV
                VA = VA + DELTVA
                N(MOLKIND) = NTRIAL
            END IF
        END IF
    END IF

    RETURN
END SUBROUTINE OUT
