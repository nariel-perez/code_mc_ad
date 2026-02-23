SUBROUTINE ESTADISTICA( NCELLMAT)
  USE SimulationData, only: CNF, RX, RY,RZ, LOCATE, N
  USE AdsorbateInput, only: NMOLEC, NATOM
  IMPLICIT NONE
  
  ! Declaración de variables
  INTEGER, INTENT(IN) :: NCELLMAT
  INTEGER :: I, IMOL, JIN, NATOMKINDI, ICNF, JCNF, KCNF
  
  ! Bucle principal
  DO I = 1, nmolec
     DO IMOL = 1, N(I)
        JIN = LOCATE(IMOL, I)
        
        DO NATOMKINDI = 1, NATOM(I)
           
           ICNF = INT(RX(JIN, NATOMKINDI, I) * NCELLMAT)
           JCNF = INT(RY(JIN, NATOMKINDI, I) * NCELLMAT)
           KCNF = INT(RZ(JIN, NATOMKINDI, I) * NCELLMAT)
           
           CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) = CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) + 1
           
        END DO
     END DO
  END DO
  
  RETURN
END SUBROUTINE ESTADISTICA

