!C------------------------------------------------------------------------------
!C    *******************************************************************
!C    ** RETURNS THE POTENTIAL ENERGY CHANGE ON ADDING AN ATOM.        **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE ADDITION **
!C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF THE ADDED ATOM   **
!C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
!C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
!C    ** REAL    SIGMA             LJ DIAMETER                         **
!C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
!C    ** REAL    RMIN              MINIMUM ALLOWED APPROACH OF ATOMS   **
!C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!C    ** DURING A TRIAL ADDITION OF AN ATOM TO THE FLUID. THE LONG     **
!C    ** RANGE CORRECTIONS IS INCLUDED.                                **
!C    *******************************************************************



SUBROUTINE POTIN(MOLKIND, SIGMA, RCUT, RMIN, DELTV, OVRLAP)

    USE SimulationData, ONLY: RX1, RY1, RZ1, NATOMKIND, NMOLEC, N, LOCATE, RX, RY, RZ, &
                             BCX, BCY, BCZ, ACELX, ACELY, ACELZ, ACEL, USS,NATOM
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: MOLKIND
    REAL :: SIGMA, RCUT, RMIN, DELTV
    LOGICAL :: OVRLAP

    ! Variables locales
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    REAL :: RXI, RYI, RZI, DELTW
    REAL :: RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
    REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
    REAL :: VLRC0, WLRC0
    INTEGER :: J, JIN
    INTEGER :: I1
    INTEGER :: IPOT, JPOT
    INTEGER :: I, K
    REAL :: RIJ
    INTEGER :: IDIST
    REAL :: XMAX, YMAX, ZMAX

    ! Constante PI
    PARAMETER (PI = 3.14159265)

    ! Inicialización de variables
    XMAX = ACELX / ACEL
    YMAX = ACELY / ACEL
    ZMAX = ACELZ / ACEL

    OVRLAP = .FALSE.
    RCUTSQ = RCUT * RCUT
    RMINSQ = RMIN * RMIN
    SIGSQ = SIGMA * SIGMA

    ! Inicialización de acumuladores
    DELTV = 0.0
    DELTW = 0.0

    ! Bucle sobre todos los átomos de la molécula
    DO I1 = 1, NATOM(MOLKIND)
        RXI = RX1(I1)
        RYI = RY1(I1)
        RZI = RZ1(I1)
        IPOT = NATOMKIND(I1, MOLKIND)

        ! Bucle sobre todas las moléculas
        DO I = 1, NMOLEC
            ! Bucle sobre todos los átomos de la molécula I
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

                    ! Verificar si hay superposición
                    IF (RIJ < RMIN) THEN
                        OVRLAP = .TRUE.
                        DELTV = 1E10
                        RETURN
                    END IF

                    ! Calcular el potencial de interacción
                    IDIST = RIJ * 1000 + 1
                    VIJ = USS(IDIST, JPOT, IPOT)

                    ! Acumular el potencial total
                    DELTV = DELTV + VIJ
                    DELTW = DELTW + WIJ
                END DO
            END DO
        END DO
    END DO

    RETURN
END SUBROUTINE POTIN
