!C---------------------------------------------------------------------------
!C---------------------------------------------------------------------------- 
!C    *******************************************************************
!C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** REAL    TEMP         TEMPERATURE                              **33
!C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
!C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
!C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
!C    ** REAL    V            POTENTIAL ENERGY                         **
!C    ** REAL    W            VIRIAL                                   **
!C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
!C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
!C    *******************************************************************
!----------------------------------------------------


SUBROUTINE CHANGE(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, GHOST, JPASOS)

    USE InputParams, ONLY: acel, acelx, acely, acelz, bcx, bcy, bcz
    USE AdsorbateInput, ONLY:  RX0, RY0, RZ0, NATOM, NMOLEC
    USE SimulationData, ONLY: RX, RY, RZ, RX1, RY1, RZ1,N, LOCATE, &
         ANX, ANGY, ANZ, &
         EXNEW, EYNEW, EZNEW
    USE RotationModule
    IMPLICIT NONE

    ! Argumentos de la subrutina
    REAL :: TEMP, Z(10), SIGMA, EPS, RCUT
    REAL :: V, VA, VG, W
    LOGICAL :: GHOST
    INTEGER :: JPASOS

    ! Variables locales
    REAL :: EXOLD, EYOLD, EZOLD
    INTEGER :: MOLKIND1, MOLKIND2
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    LOGICAL :: CREATE
    REAL :: BETA, DELTDB, DELTVA
    REAL :: RANF, DUMMY
    INTEGER :: NTRIAL, NLOC, IPULL
    REAL :: b
    REAL :: RXBE(50), RYBE(50), RZBE(50)
    REAL :: RXBE2(50), RYBE2(50), RZBE2(50)
    REAL :: RXOR, RYOR, RZOR
    REAL :: RXNEW, RYNEW, RZNEW
    INTEGER :: I
    REAL :: VOLD, VGOLD, VAOLD
    REAL :: VANT, VGANT, VAANT
    REAL :: VNUEVA, VGNUEVA, VANUEVA
    REAL :: EX, EY, EZ, RR
    REAL :: ANGUL
    PARAMETER (ANGUL = 3.14159)
    REAL :: DX, DY, DZ
    REAL :: DELTVADIN, DELTVIN
    REAL :: DELTVADOUT, DELTVOUT
    LOGICAL :: OVRLAP
    REAL :: RMIN
    REAL :: DELTCB
    REAL :: B1, B2, B3
    REAL :: rxi, ryi, rzi
    REAL :: xmax, ymax, zmax
    REAL, DIMENSION(3, 3) :: R

    ! Verificar si hay más de un tipo de molécula
    IF (NMOLEC < 2) RETURN

    ! Inicialización de variables
    xmax = ACELX / ACEL
    ymax = ACELY / ACEL
    zmax = ACELZ / ACEL

    CREATE = .FALSE.
    GHOST = .FALSE.
    BETA = 1.0 / TEMP
    RMIN = 0.75 * SIGMA

    ! Elegir tipo de molécula y recuperar posiciones
    MOLKIND1 = INT(RANF(DUMMY) * NMOLEC) + 1
    NTRIAL = N(MOLKIND1) - 1

    IF (NTRIAL < 0) RETURN

    ! Seleccionar un elemento aleatorio de la parte activa de LOCATE
    b = RANF(DUMMY)
    NLOC = INT(REAL(NTRIAL) * b) + 1
    IPULL = LOCATE(NLOC, MOLKIND1)

    ! Guardar las posiciones actuales de los átomos
    DO I = 1, NATOM(MOLKIND1)
        RXBE(I) = RX(IPULL, I, MOLKIND1)
        RYBE(I) = RY(IPULL, I, MOLKIND1)
        RZBE(I) = RZ(IPULL, I, MOLKIND1)

        RX1(I) = RXBE(I)
        RY1(I) = RYBE(I)
        RZ1(I) = RZBE(I)

        rxi = RX1(I)
        ryi = RY1(I)
        rzi = RZ1(I)
    END DO

    ! Determinar el origen
    RXOR = RXBE(1) - RX0(1, MOLKIND1)
    RYOR = RYBE(1) - RY0(1, MOLKIND1)
    RZOR = RZBE(1) - RZ0(1, MOLKIND1)

    ! Inicializar orientaciones
    EXOLD = 0
    EYOLD = 0
    EZOLD = 0
    EXNEW = 0
    EYNEW = 0
    EZNEW = 0
    IF (NATOM(MOLKIND1) > 1) THEN
        EXOLD = ANX(IPULL, MOLKIND1)
        EYOLD = ANGY(IPULL, MOLKIND1)
        EZOLD = ANZ(IPULL, MOLKIND1)
    END IF

    ! Energías viejas
    VOLD = V
    VGOLD = VG
    VAOLD = VA

    ! Calcular cambio de energía al eliminar el átomo IPULL
    CALL POTOUT(IPULL, MOLKIND1, DELTVOUT)
    CALL ADPOTOUT(IPULL, MOLKIND1, DELTVADOUT)

    ! Calcular energía antes del paso
    VANT = VOLD + DELTVOUT + DELTVADOUT
    VGANT = VGOLD + DELTVOUT
    VAANT = VAOLD + DELTVADOUT

    ! Creación de la nueva molécula
97  MOLKIND2 = INT(RANF(DUMMY) * NMOLEC) + 1
    IF (MOLKIND1 == MOLKIND2) GOTO 97

    ! Generar la posición del centro de masa del átomo de prueba
    RXNEW = RXOR
    RYNEW = RYOR
    RZNEW = RZOR

    ! Generar la nueva posición
    DO I = 1, NATOM(MOLKIND2)
        RXBE2(I) = RX0(I, MOLKIND2)
        RYBE2(I) = RY0(I, MOLKIND2)
        RZBE2(I) = RZ0(I, MOLKIND2)

        RX1(I) = RXBE2(I)
        RY1(I) = RYBE2(I)
        RZ1(I) = RZBE2(I)
    END DO

    ! Rotación molecular
    IF (NATOM(MOLKIND2) > 1) THEN
        EX = (2.0 * RANF(DUMMY) - 1.0)
        EY = (2.0 * RANF(DUMMY) - 1.0)
        EZ = (2.0 * RANF(DUMMY) - 1.0)

        RR = SQRT(EX * EX + EY * EY + EZ * EZ)
        EX = EX / RR
        EY = EY / RR
        EZ = EZ / RR

        EXNEW = EX
        EYNEW = EY
        EZNEW = EZ

        DX = ANGUL * EX
        DY = ANGUL * EY
        DZ = ANGUL * EZ

        CALL GetRotationMatrix(DX, DY, DZ, R)

        DO I = 1, NATOM(MOLKIND2)
            RX1(I) = R(1, 1) * RXBE2(I) + R(1, 2) * RYBE2(I) + R(1, 3) * RZBE2(I)
            RY1(I) = R(2, 1) * RXBE2(I) + R(2, 2) * RYBE2(I) + R(2, 3) * RZBE2(I)
            RZ1(I) = R(3, 1) * RXBE2(I) + R(3, 2) * RYBE2(I) + R(3, 3) * RZBE2(I)
        END DO
    END IF

    ! Ajustar posiciones dentro de los límites de la caja de simulación
    DO I = 1, NATOM(MOLKIND2)
        RX1(I) = RX1(I) + RXNEW
        RY1(I) = RY1(I) + RYNEW
        RZ1(I) = RZ1(I) + RZNEW

        RX1(I) = RX1(I) - BCX * xmax * ANINT(RX1(I) / xmax)
        RY1(I) = RY1(I) - BCY * ymax * ANINT(RY1(I) / ymax)
        RZ1(I) = RZ1(I) - BCZ * zmax * ANINT(RZ1(I) / zmax)

        IF (ABS(RX1(I)) > 0.5) RETURN
        IF (ABS(RY1(I)) > 0.5) RETURN
        IF (ABS(RZ1(I)) > 0.5) RETURN
    END DO

    ! Calcular el potencial de adsorción para la nueva molécula
    CALL ADPOTIN(MOLKIND2, DELTVADIN)
    IF (DELTVADIN > 100) RETURN

    ! Remover la molécula antigua
    CALL REMOVE(NLOC, IPULL, MOLKIND1)
    N(MOLKIND1) = N(MOLKIND1) - 1

    ! Calcular la nueva energía
    CALL POTIN(MOLKIND2, SIGMA, RCUT, RMIN, DELTVIN, OVRLAP)
    IF (OVRLAP) THEN
        ! Revertir cambios si hay superposición
        DO I = 1, NATOM(MOLKIND1)
            RX1(I) = RXBE(I)
            RY1(I) = RYBE(I)
            RZ1(I) = RZBE(I)
        END DO
        EXNEW = EXOLD
        EYNEW = EYOLD
        EZNEW = EZOLD
        CALL ADD(MOLKIND1)
        N(MOLKIND1) = N(MOLKIND1) + 1
        RETURN
    END IF

    ! Energía nueva
    VNUEVA = VANT + DELTVADIN + DELTVIN
    VGNUEVA = VGANT + DELTVIN
    VANUEVA = VAANT + DELTVADIN

    ! Cambio en la energía
    B1 = Z(MOLKIND2) * REAL(N(MOLKIND1) + 1)
    B2 = Z(MOLKIND1) * REAL(N(MOLKIND2) + 1)
    B3 = LOG(B1 / B2)

    DELTCB = BETA * (VNUEVA - VOLD) - B3

    ! Verificar aceptación del cambio
    IF (DELTCB < 75) THEN
        IF (DELTCB <= 0.0) THEN
            CALL ADD(MOLKIND2)
            V = VNUEVA
            VG = VGNUEVA
            VA = VANUEVA
            N(MOLKIND2) = N(MOLKIND2) + 1
        ELSE IF (EXP(-DELTCB) > RANF(DUMMY)) THEN
            CALL ADD(MOLKIND2)
            V = VNUEVA
            VG = VGNUEVA
            VA = VANUEVA
            N(MOLKIND2) = N(MOLKIND2) + 1
        ELSE
            ! Revertir cambios si no se acepta el cambio
            DO I = 1, NATOM(MOLKIND1)
                RX1(I) = RXBE(I)
                RY1(I) = RYBE(I)
                RZ1(I) = RZBE(I)
            END DO
            EXNEW = EXOLD
            EYNEW = EYOLD
            EZNEW = EZOLD
            CALL ADD(MOLKIND1)
            V = VOLD
            VG = VGOLD
            VA = VAOLD
            N(MOLKIND1) = N(MOLKIND1) + 1
        END IF
    ELSE
        ! Revertir cambios si el cambio de energía es demasiado grande
        DO I = 1, NATOM(MOLKIND1)
            RX1(I) = RXBE(I)
            RY1(I) = RYBE(I)
            RZ1(I) = RZBE(I)
        END DO
        EXNEW = EXOLD
        EYNEW = EYOLD
        EZNEW = EZOLD
        CALL ADD(MOLKIND1)
        V = VOLD
        VG = VGOLD
        VA = VAOLD
        N(MOLKIND1) = N(MOLKIND1) + 1
    END IF

    RETURN
END SUBROUTINE CHANGE
