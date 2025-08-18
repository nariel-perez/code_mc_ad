!─────────────────────────────────────────────────────────────────────
! File: Add.f90
! Agrega UNA molécula de especie MOLKIND:
!   - Busca un slot físico libre (ipull)
!   - Copia RX1/RY1/RZ1 → RX/RY/RZ en ese slot
!   - Guarda orientaciones si NATOM>1
!   - Incrementa N(MOLKIND) y actualiza LOCATE
!─────────────────────────────────────────────────────────────────────
SUBROUTINE ADD(MOLKIND)
  USE SimulationData,  ONLY : RX,RY,RZ, RX1,RY1,RZ1, ANX,ANGY,ANZ, EXNEW,EYNEW,EZNEW, &
                               LOCATE, N
  USE AdsorbateInput,  ONLY : NATOM, NMOLEC
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MOLKIND

  INTEGER :: cap, ipull, i, sp, j
  LOGICAL, ALLOCATABLE :: used(:)

  ! Capacidad máxima (nº de slots físicos disponibles)
  cap = SIZE(RX, 1)

  ! Marcar slots usados por TODAS las especies
  ALLOCATE(used(cap))
  used = .FALSE.
  DO sp = 1, NMOLEC
     IF (N(sp) > 0) THEN
        DO j = 1, N(sp)
           IF (LOCATE(j, sp) >= 1 .AND. LOCATE(j, sp) <= cap) used( LOCATE(j, sp) ) = .TRUE.
        END DO
     END IF
  END DO

  ! Buscar el primer slot libre
  ipull = 0
  DO i = 1, cap
     IF (.NOT. used(i)) THEN
        ipull = i
        EXIT
     END IF
  END DO
  DEALLOCATE(used)

  IF (ipull == 0) THEN
     WRITE(*,*) 'ADD: no free slot available (capacidad agotada).'
     STOP
  END IF

  ! Copiar la configuración trial al estado
  DO i = 1, NATOM(MOLKIND)
     RX(ipull, i, MOLKIND) = RX1(i)
     RY(ipull, i, MOLKIND) = RY1(i)
     RZ(ipull, i, MOLKIND) = RZ1(i)
  END DO

  ! Orientaciones (solo multiatómicas)
  IF (NATOM(MOLKIND) > 1) THEN
     ANX(ipull, MOLKIND) = EXNEW
     ANGY(ipull, MOLKIND) = EYNEW
     ANZ(ipull, MOLKIND) = EZNEW
  ELSE
     ANX(ipull, MOLKIND) = 0.0
     ANGY(ipull, MOLKIND) = 0.0
     ANZ(ipull, MOLKIND) = 0.0
  END IF

  ! Actualizar contadores/índices lógicos
  N(MOLKIND) = N(MOLKIND) + 1
  LOCATE(N(MOLKIND), MOLKIND) = ipull

END SUBROUTINE ADD

