MODULE EstructuraModule
  USE InputParams, only: acel, acelx, acely, acelz
  USE PhysicalConstants
  USE SimulationData
    IMPLICIT NONE
CONTAINS

    SUBROUTINE estructura(eps, nam, sigma, sigmetano, NC, diel)
        IMPLICIT NONE
        ! Variable declarations
        CHARACTER(LEN=16), INTENT(IN) :: nam
        CHARACTER(LEN=32) :: nampro
        REAL, INTENT(IN) :: eps, sigma, sigmetano, diel
        INTEGER, INTENT(INOUT) :: NC
        INTEGER :: imax, i, io_status, ipos

        ! Dynamic arrays
        REAL, ALLOCATABLE :: RXA(:), RYA(:), RZA(:)
        !INTEGER, ALLOCATABLE :: SYMBOL(:)

        ! Open output files with error handling
        OPEN(49, FILE='PRUEBA.TXT', STATUS='REPLACE', IOSTAT=io_status)
        IF (io_status /= 0) THEN
            PRINT*, "Error: No se pudo abrir PRUEBA.TXT"
            RETURN
        END IF

        OPEN(51, FILE=nam, STATUS='OLD', ACTION='READ', IOSTAT=io_status)
        IF (io_status /= 0) THEN
            PRINT*, "Error: No se pudo abrir el archivo ", TRIM(nam)
            RETURN
        END IF

        ! Read number of carbon atoms
        READ(51, *, IOSTAT=io_status) NC
        IF (io_status /= 0 .OR. NC <= 0) THEN
            PRINT*, "Error: No se pudo leer NC en ", TRIM(nam)
            CLOSE(51)
            RETURN
        END IF

        ! Allocate memory for arrays
        ALLOCATE(RXA(NC), RYA(NC), RZA(NC), SYMBOL(NC), STAT=io_status)
        IF (io_status /= 0) THEN
            PRINT*, "Error: No se pudo asignar memoria para coordenadas"
            CLOSE(51)
            RETURN
        END IF

        ! Ensure `QAC`, `RXC`, etc., are allocated
        IF (.NOT. ALLOCATED(QAC)) ALLOCATE(QAC(NC))
        IF (.NOT. ALLOCATED(SYMBOL)) ALLOCATE(SYMBOL(NC))
        IF (.NOT. ALLOCATED(RXC)) ALLOCATE(RXC(NC))
        IF (.NOT. ALLOCATED(RYC)) ALLOCATE(RYC(NC))
        IF (.NOT. ALLOCATED(RZC)) ALLOCATE(RZC(NC))
        IF (.NOT. ALLOCATED(EPSAC)) ALLOCATE(EPSAC(NC))
        IF (.NOT. ALLOCATED(SGC)) ALLOCATE(SGC(NC))

        ! Display cell dimensions
        PRINT*, "Dimensiones celda (l/2): ", ACELX/2, ACELY/2, ACELZ/2

        imax = 0
        DO i = 1, NC
            READ(51, *, IOSTAT=io_status) RXA(i), RYA(i), RZA(i), EPSAC(i), SGC(i), QAC(i), SYMBOL(i)
            IF (io_status /= 0) THEN
                PRINT*, "Error: Fallo al leer línea ", i, " en ", TRIM(nam)
                CLOSE(51)
                RETURN
            END IF

            ! Convert charge units
            QAC(i) = QAC(i) * FACTORELEC * FCLEC

            ! Check if the atom is inside the box
            IF (RXA(i) >= -ACELX/2 .AND. RXA(i) <= ACELX/2 .AND. &
                RYA(i) >= -ACELY/2 .AND. RYA(i) <= ACELY/2 .AND. &
                RZA(i) >= -ACELZ/2 .AND. RZA(i) <= ACELZ/2) THEN

                imax = imax + 1
                RXC(imax) = RXA(i) / SIGMETANO * SIGMA
                RYC(imax) = RYA(i) / SIGMETANO * SIGMA
                RZC(imax) = RZA(i) / SIGMETANO * SIGMA
                QAC(imax) = QAC(i)
                WRITE(49, *) RXC(imax), RYC(imax), RZC(imax)
            END IF
        END DO

        ! Update NC with the final count of atoms inside the box
        NC = imax
        PRINT*, "Estructura leída, ", NC, " segmentos."

        ! Close files
        CLOSE(51)
        CLOSE(49)

        ! Save truncated structure file (nombre: <base sin extensión>_truncado.txt)
        ipos = INDEX(TRIM(nam), '.', BACK=.TRUE.)
        IF (ipos > 0) THEN
            nampro = TRIM(nam(1:ipos-1)) // "_truncado.txt"
        ELSE
            nampro = TRIM(nam) // "_truncado.txt"
        END IF
        OPEN(51, FILE=nampro, STATUS='REPLACE', IOSTAT=io_status)
        IF (io_status /= 0) THEN
            PRINT*, "Error: No se pudo abrir el archivo ", TRIM(nampro)
            RETURN
        END IF

        WRITE(51, *) NC
        DO i = 1, NC
            WRITE(51, *) RXA(i), RYA(i), RZA(i), EPSAC(i), SGC(i), QAC(i) / (FCLEC * FACTORELEC), SYMBOL(i)
        END DO
        CLOSE(51)

        ! Deallocate arrays
        DEALLOCATE(RXA, RYA, RZA, SYMBOL)

    END SUBROUTINE estructura

END MODULE EstructuraModule
