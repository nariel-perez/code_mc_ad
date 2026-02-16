!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    ** SUBRUTINA POTENCIAL
!    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
!    ** UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(1000)
!------------------------------------------------------------------------
SUBROUTINE POTENCIAL(EPS, sigma, sigmetano, NC, RCUT, diel)
  USE PBC_Mod, only: rk, min_image
  USE InputParams, only: mat, cellR
  USE PhysicalConstants, only: FCLEC, FACTORELEC
  USE SimulationData
  IMPLICIT NONE
  ! Declaración de variables
  INTEGER :: NC, NKIND, INKIND, KINDI, IPOT, I, IJ, K, J
  REAL :: PI, RCELE, EPSIINKIND, SIGMINKIND, QINKIND, SIGMETANO, SIGMA, EPS
  REAL :: RXI, RYI, RZI, DELTV, SIGMA1, FACTOR, RCUTSQ, SIGSQ, SIGCUB
  REAL :: RMIN, RMINSQ, SR3, SR9, SR2, SR6
  REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, RIJ, VIJ, VIJR, ZESACT, ZI
  REAL :: diel, RCUT
  REAL(rk) :: dr(3)  ! Vector mínima imagen
  INTEGER :: IOSTAT
  
  ! Inicialización de constantes
  PI = 3.14159265
  RCELE = 0.5  ! Radio de corte potencial electrostático

  ! Abrir archivo de parámetros
  OPEN(11, FILE='LJ.dat', IOSTAT=IOSTAT)
  IF (IOSTAT /= 0) THEN
     WRITE(*, *) "Error al abrir el archivo LJ.dat"
     RETURN
  END IF
  
  ! Leer número de tipos de moléculas
    READ(11, *, IOSTAT=IOSTAT) NKIND
    IF (IOSTAT /= 0 .OR. NKIND <= 0) THEN
       WRITE(*, *) "Error al leer NKIND"
       CLOSE(11)
       RETURN
    END IF

    ! Leer parámetros para cada tipo de molécula
    DO INKIND = 1, NKIND
       READ(11, *, IOSTAT=IOSTAT) KINDI, EPSIINKIND, SIGMINKIND, QINKIND
       IF (IOSTAT /= 0) THEN
          WRITE(*, *) "Error al leer parámetros para el tipo de molécula", INKIND
          CLOSE(11)
          RETURN
       END IF
       
        EPSI(INKIND) = EPSIINKIND
        SIGM(INKIND) = SIGMINKIND
        Q(INKIND) = QINKIND
        Q(INKIND) = Q(INKIND) * FACTORELEC
        Q(INKIND) = Q(INKIND) * FCLEC
     END DO
     CLOSE(11)
     
     ! Calcular potenciales
     DO IPOT = 1, NKIND
        WRITE(*,*) 'POTENCIAL PARA ', IPOT,' DE ', NKIND
        DO I = -mat/2, mat/2
           RZI = REAL(I) / REAL(mat)
           DO IJ = -mat/2, mat/2
              RYI = REAL(IJ) / REAL(mat)
              DO K = -mat/2, mat/2
                 RXI = REAL(K) / REAL(mat)
                 DELTV = 0.0
                 VIJ = 0.0
                 VIJR = 0.0
                 
                 DO J = 1, NC
                    SIGMA1 = SIGMA * (SGC(J) + SIGM(IPOT)) / (2 * SIGMETANO)
                    FACTOR = SQRT(EPSI(IPOT) * EPSAC(J)) / EPS
                    RCUTSQ = RCUT * RCUT
                    SIGSQ = SIGMA1 * SIGMA1
                    SIGCUB = SIGSQ * SIGMA1
                    RMIN = 0.5 * SIGMA1
                    RMINSQ = RMIN * RMIN
                    SR3 = (SIGMA1 / RCUT) ** 3
                    SR9 = SR3 ** 3
                    
                    ! Vector mínima imagen usando cellR (unidades reducidas)
                    dr = min_image(cellR, [real(RXC(J),rk), real(RYC(J),rk), real(RZC(J),rk)], &
                                          [real(RXI,rk), real(RYI,rk), real(RZI,rk)])
                    RXIJ = real(dr(1))
                    RYIJ = real(dr(2))
                    RZIJ = real(dr(3))
                    RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                    
                    IF (RIJSQ < RMINSQ) THEN
                       VIJ = 1E6
                       DELTV = 1E6
                       EXIT
                    ELSE
                       SR2 = SIGSQ / RIJSQ
                       SR6 = SR2 * SR2 * SR2
                       VIJ = SR6 * (SR6 - 1.0) * FACTOR
                       IF (RIJSQ > (8 * RMIN)**2) VIJ = 0.0
                       
                       IF (RIJSQ < RCELE**2) THEN
                          RIJ = SQRT(RIJSQ)
                          ZESACT = QAC(J)
                          ZI = Q(IPOT)
                          VIJR = ZI * ZESACT * (1 / RIJ - 1 / RCELE + (1 / RCELE**2) * (RIJ - RCELE))
                       ELSE
                          VIJR = 0.0
                       END IF
                    END IF
                    
                    DELTV = DELTV + 4 * VIJ + VIJR
                    END DO
                    
                    UADS(K, IJ, I, IPOT) = DELTV
                    IF (DELTV > 0) DELTV = 0.0
                    DELTV = DELTV / EPS * 8.31
                 END DO
              END DO
              !WRITE(*,*) DELTV,K,IJ,I,IPOT 
           END DO
           WRITE(*,*) 'FIN DE POTENCIAL PARA ', IPOT,' DE ', NKIND
        END DO
        
        WRITE(*, *) 'Energia Calculada'
        RETURN
      END SUBROUTINE POTENCIAL
      
