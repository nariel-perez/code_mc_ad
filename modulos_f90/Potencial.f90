!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    ** SUBRUTINA POTENCIAL
!    ** CALCULA EL POTENCIAL EN CADA PUNTO Y LO GUARDA EN LA MATRIZ
!    ** UADS(-225:225,-225:225,-225:225),USS(1000),FLAG(1000)
!------------------------------------------------------------------------
SUBROUTINE POTENCIAL(EPS, sigma, sigmetano, NC, RCUT, diel)
  USE InputParams, only: acel,acelx, acely, acelz,bcx, bcy,bcz,mat
  USE PhysicalConstants, only: FCLEC, FACTORELEC
  USE SimulationData
  IMPLICIT NONE
  
  ! Declaración de variables
  INTEGER :: MOLKIND, NC, NKIND, INKIND, KINDI, IPOT, I, IJ, K, J
  REAL :: PI, RCELE, EPSIINKIND, SIGMINKIND, QINKIND, SIGMETANO, SIGMA, EPS
  REAL :: RXI, RYI, RZI, DELTV, SIGMA1, FACTOR, RCUTSQ, SIGSQ, SIGCUB
  REAL :: RMIN, RMINSQ, SR3, SR9, SR2, SR6, VLRC0, WLRC0, DELTW
  REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, RIJ, VIJ, WIJ, VIJR, ZESACT, ZI
  REAL :: diel, xmax, ymax, zmax, REDELEC, RCUT
  INTEGER :: IOSTAT
  
  ! Inicialización de constantes
  PI = 3.14159265
  RCELE = 0.5  ! Radio de corte potencial electrostático
  
  ! Calcular dimensiones máximas
  xmax = acelx / acel
  ymax = acely / acel
  zmax = acelz / acel

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
                    
                    RXIJ = RXI - RXC(J)
                    RYIJ = RYI - RYC(J)
                    RZIJ = RZI - RZC(J)
                    
                    RXIJ = RXIJ - BCX * xmax * ANINT(RXIJ / xmax)
                    RYIJ = RYIJ - BCY * ymax * ANINT(RYIJ / ymax)
                    RZIJ = RZIJ - BCZ * zmax * ANINT(RZIJ / zmax)
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
                       
                       WIJ = SR6 * (SR6 - 0.5) * FACTOR
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
                    DELTW = DELTW + WIJ
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
      
