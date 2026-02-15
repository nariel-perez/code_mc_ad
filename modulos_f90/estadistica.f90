SUBROUTINE ESTADISTICA( NCELLMAT)
  USE PBC_Mod, only: rk, cart_to_frac
  USE InputParams, only: cellR
  USE SimulationData, only: CNF, RX, RY,RZ, LOCATE, N
  USE AdsorbateInput, only: NMOLEC, NATOM
  IMPLICIT NONE

  ! Declaración de variables
  INTEGER, INTENT(IN) :: NCELLMAT
  INTEGER :: I, IMOL, JIN, NATOMKINDI, ICNF, JCNF, KCNF, half
  real(rk) :: pos(3), s(3)

  half = NCELLMAT / 2

  ! Bucle principal
  DO I = 1, nmolec
     DO IMOL = 1, N(I)
        JIN = LOCATE(IMOL, I)

        DO NATOMKINDI = 1, NATOM(I)

           pos = [real(RX(JIN, NATOMKINDI, I), rk), &
                  real(RY(JIN, NATOMKINDI, I), rk), &
                  real(RZ(JIN, NATOMKINDI, I), rk)]
           s = cart_to_frac(cellR, pos)

           ICNF = INT(s(1) * NCELLMAT)
           JCNF = INT(s(2) * NCELLMAT)
           KCNF = INT(s(3) * NCELLMAT)

           ! Clamp to valid range
           ICNF = max(-half, min(half, ICNF))
           JCNF = max(-half, min(half, JCNF))
           KCNF = max(-half, min(half, KCNF))

           CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) = CNF(ICNF, JCNF, KCNF, I, NATOMKINDI) + 1

        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE ESTADISTICA

