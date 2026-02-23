!C-----------------------------------------------------------------------------
!C-----------------------------------------------------------------------------
!C    *******************************************************************
!C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
!C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
!C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
!C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
!C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
!C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
!C    ** REAL    SIGMA             LJ DIAMETER                         **
!C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
!C    ** RANGE CORRECTIONS IS INCLUDED.                                **
!C    *******************************************************************
!----------------------------------------------------


SUBROUTINE POTOUT(IPULL, MOLKIND, DELTV)
    USE InputParams, ONLY:  acel, acelx, acely, acelz, bcx, bcy, bcz
    USE AdsorbateInput, ONLY: NATOMKIND, NMOLEC,NATOM
    USE SimulationData, ONLY: RX, RY, RZ, N, LOCATE, &
                              USS
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: IPULL, MOLKIND
    REAL :: DELTV

    ! Variables locales
    REAL :: RXI, RYI, RZI
    REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, DELTW
    INTEGER :: J, JIN
    INTEGER :: I1
    INTEGER :: IPOT, JPOT
    INTEGER :: I, K
    REAL :: RIJ
    INTEGER :: IDIST
    REAL :: XMAX, YMAX, ZMAX

    ! Inicialización de variables
    XMAX = ACELX / ACEL
    YMAX = ACELY / ACEL
    ZMAX = ACELZ / ACEL

    ! Inicialización de acumuladores
    DELTV = 0.0
    DELTW = 0.0

    ! Bucle sobre todos los átomos de la molécula
    DO I1 = 1, NATOM(MOLKIND)
        RXI = RX(IPULL, I1, MOLKIND)
        RYI = RY(IPULL, I1, MOLKIND)
        RZI = RZ(IPULL, I1, MOLKIND)
        IPOT = NATOMKIND(I1, MOLKIND)

        ! Bucle sobre todas las moléculas
        DO I = 1, NMOLEC
            IF (I /= MOLKIND) THEN
                ! Bucle sobre todas las moléculas de tipo I
                DO J = 1, N(I)
                    ! Seleccionar átomos activos del array LOCATE
                    JIN = LOCATE(J, I)

                    ! Bucle sobre todos los átomos de la molécula I
                    DO K = 1, NATOM(I)
                        JPOT = NATOMKIND(K, I)

                        ! Calcular distancias entre átomos
                        RXIJ = RXI - RX(JIN, K, I)
                        RYIJ = RYI - RY(JIN, K, I)
                        RZIJ = RZI - RZ(JIN, K, I)

                        ! Ajustar distancias dentro de los límites de la caja de simulación
                        RXIJ = RXIJ - BCX * XMAX * ANINT(RXIJ / XMAX)
                        RYIJ = RYIJ - BCY * YMAX * ANINT(RYIJ / YMAX)
                        RZIJ = RZIJ - BCZ * ZMAX * ANINT(RZIJ / ZMAX)

                        ! Calcular la distancia al cuadrado
                        RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                        RIJ = SQRT(RIJSQ)

                        ! Calcular el potencial de interacción
                        IDIST = INT(RIJ * 1000.0 + 1.0)
                        VIJ = USS(IDIST, JPOT, IPOT)

                        ! Acumular el potencial total
                        DELTV = DELTV + VIJ
                        DELTW = DELTW + WIJ
                    END DO
                END DO
            ELSE
                ! Bucle sobre todas las moléculas del mismo tipo (MOLKIND), porque fueron excluidad en el bucle anterior
                DO J = 1, N(I)
                    ! Seleccionar átomos activos del array LOCATE
                    JIN = LOCATE(J, I)

                    ! Evitar el átomo que se está eliminando (IPULL)
                    IF (JIN /= IPULL) THEN
                        ! Bucle sobre todos los átomos de la molécula I
                        DO K = 1, NATOM(I)
                            JPOT = NATOMKIND(K, I)

                            ! Calcular distancias entre átomos
                            RXIJ = RXI - RX(JIN, K, I)
                            RYIJ = RYI - RY(JIN, K, I)
                            RZIJ = RZI - RZ(JIN, K, I)

                            ! Ajustar distancias dentro de los límites de la caja de simulación
                            RXIJ = RXIJ - BCX * XMAX * ANINT(RXIJ / XMAX)
                            RYIJ = RYIJ - BCY * YMAX * ANINT(RYIJ / YMAX)
                            RZIJ = RZIJ - BCZ * ZMAX * ANINT(RZIJ / ZMAX)

                            ! Calcular la distancia al cuadrado
                            RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                            RIJ = SQRT(RIJSQ)

                            ! Calcular el potencial de interacción
                            IDIST = INT(RIJ * 1000.0 + 1.0)
                            VIJ = USS(IDIST, JPOT, IPOT)

                            ! Acumular el potencial total
                            DELTV = DELTV + VIJ
                            DELTW = DELTW + WIJ
                        END DO
                    END IF
                END DO
            END IF
        END DO
    END DO

    ! Cambiar el signo de DELTV y DELTW para una eliminación
    DELTV = -DELTV
    DELTW = -DELTW

    RETURN
END SUBROUTINE POTOUT
