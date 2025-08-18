!─────────────────────────────────────────────────────────────────────
! File: Remove.f90
! Elimina la molécula #NLOC de la especie MOLKIND:
!   - ipull es el slot físico a liberar
!   - Compacta LOCATE a partir de NLOC
!   - Decrementa N(MOLKIND)
!   - Limpia el último LOCATE y deja libre ipull
!─────────────────────────────────────────────────────────────────────
SUBROUTINE REMOVE(NLOC, IPULL, MOLKIND)
  USE SimulationData, ONLY : LOCATE, N, ANX,ANGY,ANZ
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NLOC, IPULL, MOLKIND
  INTEGER :: k

  IF (N(MOLKIND) <= 0) RETURN
  IF (NLOC < 1 .OR. NLOC > N(MOLKIND)) RETURN

  ! Compactar LOCATE: correr hacia arriba desde NLOC
  DO k = NLOC + 1, N(MOLKIND)
     LOCATE(k-1, MOLKIND) = LOCATE(k, MOLKIND)
  END DO

  ! Decrementar cantidad lógica y limpiar el “último”
  LOCATE(N(MOLKIND), MOLKIND) = 0
  N(MOLKIND) = N(MOLKIND) - 1

  ! (opcional) limpiar orientaciones del slot físico liberado
  ANX(IPULL, MOLKIND) = 0.0
  ANGY(IPULL, MOLKIND) = 0.0
  ANZ(IPULL, MOLKIND) = 0.0
END SUBROUTINE REMOVE

