!C--------------------------------------------------------------------
!C    *******************************************************************
!C    ** ROUTINE TO ATTEMPT A TRIAL CREATION                           **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** REAL    TEMP            TEMPERATURE                           **
!C    ** REAL    Z               ABSOLUTE ACTIVITY                     **
!C    ** REAL    SIGMA           LENNARD-JONES DIAMETER                **
!C    ** REAL    RCUT            CUT-OFF DISTANCE                      **
!C    ** REAL    V               POTENTIAL ENERGY                      **
!C    ** REAL    W               VIRIAL                                **
!C    ** INTEGER N               NUMBER OF ATOMS BEFORE TRIAL CREATION **
!C    ** LOGICAL CREATE          TRUE FOR A SUCCESSFUL CREATION        **
!C    *******************************************************************
      
SUBROUTINE IN(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, CREATE, CR, JPASOS, CANONICALMOLECULES)

  USE InputParams, ONLY: acel, acelx, acely, acelz, bcx, bcy, bcz
  USE AdsorbateInput, ONLY: RX0, RY0, RZ0, NMOLEC,NATOM
  USE SimulationData, ONLY: RX1, RY1, RZ1, EXNEW, EYNEW, EZNEW, N
  USE RotationModule
  
  IMPLICIT NONE
  
  ! Argumentos de la subrutina
  INTEGER :: JPASOS, CANONICALMOLECULES
  REAL :: TEMP, Z(10), SIGMA, EPS, RCUT
  REAL :: V, VA, VG, W, CR
  LOGICAL :: CREATE
  
  ! Variables locales
  INTEGER :: NMAX
  PARAMETER (NMAX = 15000)
  REAL :: ANGUL
  PARAMETER (ANGUL = 3.14159)
  REAL :: RXBE(50), RYBE(50), RZBE(50)
  REAL :: BETA, RXNEW, RYNEW, RZNEW
  REAL :: RANF, DUMMY, RMIN
  INTEGER :: NTRIAL
  LOGICAL :: OVRLAP
  REAL :: XMAX, YMAX, ZMAX
  INTEGER :: MOLKIND
  INTEGER :: I
  REAL :: DX, DY, DZ
  REAL :: EX, EY, EZ, RR
  REAL :: DELTVA, DELTV, DELTW, DELTCB
  REAL, DIMENSION(3, 3) :: R
  
  ! Inicialización de variables
  EXNEW = 0.0
  EYNEW = 0.0
  EZNEW = 0.0
  
  XMAX = ACELX / ACEL
  YMAX = ACELY / ACEL
  ZMAX = ACELZ / ACEL
  
  CREATE = .FALSE.
  BETA = 1.0 / TEMP
  RMIN = 0.75 * SIGMA
  DELTW = 0.0
  
  ! Elegir tipo de molécula
  MOLKIND = INT(RANF(DUMMY) * NMOLEC) + 1
  NTRIAL = N(MOLKIND) + 1
  
  IF (NTRIAL >= NMAX) RETURN
  
  ! Generar la posición del centro de masa del átomo de prueba
  RXNEW = (RANF(DUMMY) - 0.5) * XMAX
  RYNEW = (RANF(DUMMY) - 0.5) * YMAX
  RZNEW = (RANF(DUMMY) - 0.5) * ZMAX
  
  ! Generar la nueva posición
  DO I = 1, NATOM(MOLKIND)
     RXBE(I) = RX0(I, MOLKIND)
     RYBE(I) = RY0(I, MOLKIND)
     RZBE(I) = RZ0(I, MOLKIND)
     RX1(I) = RXBE(I)
     RY1(I) = RYBE(I)
     RZ1(I) = RZBE(I)
  END DO
  
  ! Rotación molecular
  IF (NATOM(MOLKIND) > 1) THEN
     EX = (2.0 * RANF(DUMMY) - 1.0)
     EY = (2.0 * RANF(DUMMY) - 1.0)
     EZ = (2.0 * RANF(DUMMY) - 1.0)

     RR = SQRT(EX**2 + EY**2 + EZ**2)
     EX = EX / RR
     EY = EY / RR
     EZ = EZ / RR
     
     DX = ANGUL * EX
     DY = ANGUL * EY
     DZ = ANGUL * EZ
     
     CALL GetRotationMatrix(DX, DY, DZ, R)
     
     DO I = 1, NATOM(MOLKIND)
        RX1(I) = R(1, 1) * RXBE(I) + R(1, 2) * RYBE(I) + R(1, 3) * RZBE(I)
        RY1(I) = R(2, 1) * RXBE(I) + R(2, 2) * RYBE(I) + R(2, 3) * RZBE(I)
        RZ1(I) = R(3, 1) * RXBE(I) + R(3, 2) * RYBE(I) + R(3, 3) * RZBE(I)
     END DO
  END IF
      
  ! Ajustar posiciones dentro de los límites de la caja de simulación
  DO I = 1, NATOM(MOLKIND)
     RX1(I) = RX1(I) + RXNEW
     RY1(I) = RY1(I) + RYNEW
     RZ1(I) = RZ1(I) + RZNEW
     
     RX1(I) = RX1(I) - BCX * XMAX * ANINT(RX1(I) / XMAX)
     RY1(I) = RY1(I) - BCY * YMAX * ANINT(RY1(I) / YMAX)
     RZ1(I) = RZ1(I) - BCZ * ZMAX * ANINT(RZ1(I) / ZMAX)
     
     IF (ABS(RX1(I)) > 0.5) RETURN
     IF (ABS(RY1(I)) > 0.5) RETURN
     IF (ABS(RZ1(I)) > 0.5) RETURN
  END DO
  
  ! Calcular el potencial de adsorción
  CALL ADPOTIN(MOLKIND, DELTVA)
  IF (DELTVA > 1000) RETURN
  
  ! Calcular el potencial de interacción
  CALL POTIN(MOLKIND, SIGMA, RCUT, RMIN, DELTV, OVRLAP)
  IF (OVRLAP) RETURN
  
  ! Verificar la aceptación del cambio
  IF (.NOT. OVRLAP) THEN
     DELTCB = BETA * (DELTV + DELTVA) - LOG(Z(MOLKIND) / REAL(NTRIAL))
     
     IF (DELTCB < 75.0) THEN
        IF (DELTCB <= 0.0) THEN
           CREATE = .TRUE.
           CALL ADD(MOLKIND)
           V = V + DELTV + DELTVA
           VG = VG + DELTV
           VA = VA + DELTVA
           W = W + DELTW
           N(MOLKIND) = NTRIAL
        ELSE IF (EXP(-DELTCB) > RANF(DUMMY)) THEN
           CREATE = .TRUE.
           CALL ADD(MOLKIND)
           CR = CR + 1
           V = V + DELTV + DELTVA
           VG = VG + DELTV
           VA = VA + DELTVA
           W = W + DELTW
           N(MOLKIND) = NTRIAL
        END IF
     END IF
  END IF
  
  RETURN
END SUBROUTINE IN

