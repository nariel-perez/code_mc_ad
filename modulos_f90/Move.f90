!C-------------------------------------------------------------------------
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


SUBROUTINE MOVE(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, GHOST, JPASOS)

    USE PBC_Mod, ONLY: rk, cart_to_frac, frac_to_cart
    USE GeomUtils, ONLY: wrap_by_pbc
    USE InputParams, ONLY: cellR
    USE AdsorbateInput, ONLY: RX0, RY0, RZ0, NATOM, NMOLEC
    USE SimulationData, ONLY: RX, RY, RZ, RX1, RY1, RZ1,N, LOCATE, ANX, ANGY, ANZ, &
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
    REAL :: DELTVADOUT2, DELTVOUT2
    INTEGER :: MOLKIND
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    LOGICAL :: CREATE
    REAL :: BETA, DELTDB, DELTVA
    REAL :: RANF, DUMMY
    INTEGER :: NTRIAL, NLOC, IPULL
    REAL :: b
    REAL :: DELTAX, DELTAY, DELTAZ, DELTA
    REAL :: RXBE(50), RYBE(50), RZBE(50)
    REAL :: RXBE2(50), RYBE2(50), RZBE2(50)
    INTEGER :: I, IJ
    REAL :: VOLD, VGOLD, VAOLD
    REAL :: VANT, VGANT, VAANT
    REAL :: VNUEVA, VGNUEVA, VANUEVA
    REAL :: EX, EY, EZ, RR
    REAL :: DX, DY, DZ
    REAL :: DELTVADOUT, DELTVOUT
    REAL :: RMIN
    REAL :: DELTCB
    REAL :: rxnew, rynew, rznew
    REAL :: rxi, ryi, rzi
    REAL, DIMENSION(3, 3) :: R
    REAL(rk) :: pos(3), s(3)  ! Para wrap de posiciones

    CREATE = .FALSE.
    GHOST = .FALSE.
    BETA = 1.0 / TEMP
    RMIN = 0.75 * SIGMA

    ! Elegir tipo de molécula y recuperar posiciones
    MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1
    NTRIAL = N(MOLKIND) - 1

    IF (NTRIAL < 0) RETURN

    ! Seleccionar un elemento aleatorio de la parte activa de LOCATE
    b = RANF(DUMMY)
    NLOC = INT(REAL(NTRIAL) * b) + 1
    IPULL = LOCATE(NLOC, MOLKIND)

    ! Guardar las posiciones actuales de los átomos
    DO I = 1, NATOM(MOLKIND)
        RXBE(I) = RX(IPULL, I, MOLKIND)
        RYBE(I) = RY(IPULL, I, MOLKIND)
        RZBE(I) = RZ(IPULL, I, MOLKIND)

        RX1(I) = RXBE(I)
        RY1(I) = RYBE(I)
        RZ1(I) = RZBE(I)

        rxi = RX1(I)
        ryi = RY1(I)
        rzi = RZ1(I)
    END DO

    ! Inicializar orientaciones
    EXOLD = 0
    EYOLD = 0
    EZOLD = 0
    EXNEW = 0
    EYNEW = 0
    EZNEW = 0
    IF (NATOM(MOLKIND) > 1) THEN
        EXOLD = ANX(IPULL, MOLKIND)
        EYOLD = ANGY(IPULL, MOLKIND)
        EZOLD = ANZ(IPULL, MOLKIND)
    END IF

    ! Energías viejas
    VOLD = V
    VGOLD = VG
    VAOLD = VA

    ! Calcular cambio de energía al eliminar el átomo IPULL
    CALL POTOUT(IPULL, MOLKIND, DELTVOUT)
    CALL ADPOTOUT(IPULL, MOLKIND, DELTVADOUT)

    ! Calcular energía antes del paso
    VANT = VOLD + DELTVOUT + DELTVADOUT
    VGANT = VGOLD + DELTVOUT
    VAANT = VAOLD + DELTVADOUT

    ! Elegir entre rotación o traslación
    IJ = INT(RANF(DUMMY) * 2) + 1
    IF (NATOM(MOLKIND) == 1) IJ = 2

    SELECT CASE (IJ)
        CASE (1)
            ! Rotación molecular
            IF (NATOM(MOLKIND) > 1) THEN
                DO I = 1, NATOM(MOLKIND)
                    RXBE2(I) = RX0(I, MOLKIND)
                    RYBE2(I) = RY0(I, MOLKIND)
                    RZBE2(I) = RZ0(I, MOLKIND)
                END DO

                rxnew = RX1(2)
                rynew = RY1(2)
                rznew = RZ1(2)

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

                DX = 3.14159 * EX
                DY = 3.14159 * EY
                DZ = 3.14159 * EZ

                CALL GetRotationMatrix(DX, DY, DZ, R)

                DO I = 1, NATOM(MOLKIND)
                    RX1(I) = R(1, 1) * RXBE2(I) + R(1, 2) * RYBE2(I) + R(1, 3) * RZBE2(I)
                    RY1(I) = R(2, 1) * RXBE2(I) + R(2, 2) * RYBE2(I) + R(2, 3) * RZBE2(I)
                    RZ1(I) = R(3, 1) * RXBE2(I) + R(3, 2) * RYBE2(I) + R(3, 3) * RZBE2(I)
                END DO

                DO I = 1, NATOM(MOLKIND)
                    RX1(I) = RX1(I) + rxnew
                    RY1(I) = RY1(I) + rynew
                    RZ1(I) = RZ1(I) + rznew

                    ! Wrap usando cellR (respeta flags PBC)
                    pos = [real(RX1(I),rk), real(RY1(I),rk), real(RZ1(I),rk)]
                    s = cart_to_frac(cellR, pos)
                    call wrap_by_pbc(s(1), s(2), s(3), cellR%pbc(1), cellR%pbc(2), cellR%pbc(3))
                    pos = frac_to_cart(cellR, s)
                    RX1(I) = real(pos(1))
                    RY1(I) = real(pos(2))
                    RZ1(I) = real(pos(3))

                    IF (ABS(RX1(I)) > 0.5) RETURN
                    IF (ABS(RY1(I)) > 0.5) RETURN
                    IF (ABS(RZ1(I)) > 0.5) RETURN
                END DO
            END IF

        CASE (2)
            ! Traslación
            DELTA = 0.01
            EXNEW = EXOLD
            EYNEW = EYOLD
            EZNEW = EZOLD

            DELTAX = (RANF(DUMMY) - 0.5) * SIGMA * DELTA
            DELTAY = (RANF(DUMMY) - 0.5) * SIGMA * DELTA
            DELTAZ = (RANF(DUMMY) - 0.5) * SIGMA * DELTA

            DO I = 1, NATOM(MOLKIND)
                RX1(I) = RX1(I) + DELTAX
                RY1(I) = RY1(I) + DELTAY
                RZ1(I) = RZ1(I) + DELTAZ

                ! Wrap usando cellR (respeta flags PBC)
                pos = [real(RX1(I),rk), real(RY1(I),rk), real(RZ1(I),rk)]
                s = cart_to_frac(cellR, pos)
                call wrap_by_pbc(s(1), s(2), s(3), cellR%pbc(1), cellR%pbc(2), cellR%pbc(3))
                pos = frac_to_cart(cellR, s)
                RX1(I) = real(pos(1))
                RY1(I) = real(pos(2))
                RZ1(I) = real(pos(3))

                IF (ABS(RX1(I)) > 0.5) RETURN
                IF (ABS(RY1(I)) > 0.5) RETURN
                IF (ABS(RZ1(I)) > 0.5) RETURN
            END DO
    END SELECT

    ! Actualizar posiciones y orientaciones
    DO I = 1, NATOM(MOLKIND)
        RX(IPULL, I, MOLKIND) = RX1(I)
        RY(IPULL, I, MOLKIND) = RY1(I)
        RZ(IPULL, I, MOLKIND) = RZ1(I)
    END DO
    ANX(IPULL, MOLKIND) = EXNEW
    ANGY(IPULL, MOLKIND) = EYNEW
    ANZ(IPULL, MOLKIND) = EZNEW

    ! Calcular cambio de energía después del movimiento
    CALL POTOUT(IPULL, MOLKIND, DELTVOUT2)
    CALL ADPOTOUT(IPULL, MOLKIND, DELTVADOUT2)

    ! Energía nueva
    VNUEVA = VANT - DELTVOUT2 - DELTVADOUT2
    VGNUEVA = VGANT - DELTVOUT2
    VANUEVA = VAANT - DELTVADOUT2

    ! Cambio en la energía
    DELTCB = BETA * (VNUEVA - VOLD)

    ! Verificar aceptación del movimiento
    IF (DELTCB < 75) THEN
        IF (DELTCB <= 0.0) THEN
            V = VNUEVA
            VG = VGNUEVA
            VA = VANUEVA
        ELSE IF (EXP(-DELTCB) > RANF(DUMMY)) THEN
            V = VNUEVA
            VG = VGNUEVA
            VA = VANUEVA
        ELSE
            ! Revertir cambios si no se acepta el movimiento
            DO I = 1, NATOM(MOLKIND)
                RX(IPULL, I, MOLKIND) = RXBE(I)
                RY(IPULL, I, MOLKIND) = RYBE(I)
                RZ(IPULL, I, MOLKIND) = RZBE(I)
            END DO
            ANX(IPULL, MOLKIND) = EXOLD
            ANGY(IPULL, MOLKIND) = EYOLD
            ANZ(IPULL, MOLKIND) = EZOLD

            V = VOLD
            VG = VGOLD
            VA = VAOLD
        END IF
    ELSE
        ! Revertir cambios si el cambio de energía es demasiado grande
        DO I = 1, NATOM(MOLKIND)
            RX(IPULL, I, MOLKIND) = RXBE(I)
            RY(IPULL, I, MOLKIND) = RYBE(I)
            RZ(IPULL, I, MOLKIND) = RZBE(I)
        END DO
        ANX(IPULL, MOLKIND) = EXOLD
        ANGY(IPULL, MOLKIND) = EYOLD
        ANZ(IPULL, MOLKIND) = EZOLD

        V = VOLD
        VG = VGOLD
        VA = VAOLD
    END IF

    RETURN
END SUBROUTINE MOVE
