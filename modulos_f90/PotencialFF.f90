!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    ** SUBRUTINA POTENCIAL2
!    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
!    ** USS(1000),FLAG(1000)
!------------------------------------------------------------------------

SUBROUTINE POTENCIALFF(EPS, sigma, sigmetano, NC, RCUT, diel)
    USE PhysicalConstants, only: FCLEC, FACTORELEC
    USE SimulationData
    IMPLICIT NONE

    ! Declaración de variables
    INTEGER :: MOLKIND, NKIND, INKIND, KINDI, IPOT, JPOT, I, NC
    REAL :: PI, RCELE, SIGMA1, FACTOR, RCUTSQ, SIGSQ, SIGCUB, RMIN, RMINSQ
    REAL :: SR3, SR9, SR2, SR6, VLRC0, WLRC0, DELTV, DELTW, RZI, RIJSQ, VIJ, WIJ, VIJE1
    REAL :: EPS, sigma, sigmetano, RCUT, diel, REDELEC
    INTEGER :: IOSTAT
    ! Inicialización de constantes
    PI = 3.14159265
    RCELE = 0.35  ! Radio de corte potencial electrostático

    
    ! Abrir archivo de parámetros
    OPEN(11, FILE='LJ.dat', IOSTAT=IOSTAT)
    IF (IOSTAT /= 0) THEN
        WRITE(*, *) "Error al abrir el archivo LJ.dat"
        RETURN
    END IF

    ! Leer número de tipos de moléculas y factor de reducción electrostática
    READ(11, *, IOSTAT=IOSTAT) NKIND, REDELEC
    IF (IOSTAT /= 0 .OR. NKIND <= 0) THEN
        WRITE(*, *) "Error al leer NKIND o REDELEC"
        CLOSE(11)
        RETURN
    END IF

    ! Escribir encabezados
    WRITE(*, '(A10, A10, A10, A10)') 'TYPE', 'EPSILON', 'SIGMA', 'CHARGE'
    WRITE(*, '(A10, A10, A10, A10)') '----', '-------', '-----', '------'

    ! Leer parámetros para cada tipo de molécula
    DO INKIND = 1, NKIND
        READ(11, *, IOSTAT=IOSTAT) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)
        IF (IOSTAT /= 0) THEN
            WRITE(*, *) "Error al leer parámetros para el tipo de molécula", INKIND
            CLOSE(11)
            RETURN
        END IF
        WRITE(*, *) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)

        ! Ajustar carga electrostática
        Q(INKIND) = Q(INKIND) * FACTORELEC
        Q(INKIND) = REAL(Q(INKIND)) * FCLEC
    END DO
    CLOSE(11)

    ! Calcular potenciales

    DO IPOT = 1, NKIND
        DO JPOT = IPOT, NKIND  ! Solo calcular para JPOT >= IPOT
            SIGMA1 = SIGMA * (SIGM(IPOT) + SIGM(JPOT)) / (2 * SIGMETANO)
            FACTOR = SQRT(EPSI(IPOT) * EPSI(JPOT)) / EPS
            RCUTSQ = RCUT * RCUT
            SIGSQ = SIGMA1 * SIGMA1
            SIGCUB = SIGSQ * SIGMA1
            RMIN = 0.5 * SIGMA1
            RMINSQ = RMIN * RMIN

            ! Corrección de largo alcance para Lennard-Jones
            SR3 = (SIGMA1 / RCUT) ** 3
            SR9 = SR3 ** 3
            VLRC0 = (8.0 / 9.0) * PI * SIGCUB * (SR9 - 3.0 * SR3)
            WLRC0 = (16.0 / 9.0) * PI * SIGCUB * (2.0 * SR9 - 3.0 * SR3)

            ! Inicializar acumuladores
            DELTV = 0.0
            DELTW = 0.0

            ! Calcular potencial en cada punto
            DO I = 1, 5000
                RZI = REAL(I) / 1000.0
                RIJSQ = RZI * RZI

                IF (RIJSQ < RMINSQ) THEN
                    VIJ = 1E10 * SIGM(IPOT) * SIGM(JPOT)
                ELSE
                    SR2 = SIGSQ / RIJSQ
                    SR6 = SR2 * SR2 * SR2
                    VIJ = SR6 * (SR6 - 1.0)
                    IF (RZI > 8 * RMIN) VIJ = 0
                    WIJ = SR6 * (SR6 - 0.5)
                END IF

                DELTV = 4.0 * VIJ * FACTOR

                ! Potencial electrostático
                VIJE1 = Q(IPOT) * Q(JPOT) * (1 / RZI - 1 / RCELE + (1 / RCELE**2) * (RZI - RCELE))
                IF (RZI > RCELE) VIJE1 = 0
                IF (RZI < RMIN) VIJE1 = 1E9 * ABS(Q(IPOT) * Q(JPOT))

                ! Almacenar en la matriz USS (simetría aplicada)
                USS(I, IPOT, JPOT) = DELTV + VIJE1
                USS(I, JPOT, IPOT) = USS(I, IPOT, JPOT)  ! Asignación simétrica
            END DO
        END DO
    END DO

    ! Cerrar archivo y liberar memoria
    WRITE(*, *) '-----------------------'

    RETURN
END SUBROUTINE POTENCIALFF
