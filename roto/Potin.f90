!─────────────────────────────────────────────────────────────────────
! File: Potin.f90         (Fortran 2003+)
! Propósito:
!   Energía de interacción de la molécula de PRUEBA (RX1/RY1/RZ1,
!   especie MOLKIND) con TODAS las moléculas existentes, usando USS.
!   Versión general para celdas triclinic y PBC parciales.
!
! Notas:
!   - Distancias en fraccionales de la celda "patrón" (cellR = A/ACEL).
!   - r = sqrt( s^T G s ), con G = A^T A (métrica) precomputada localmente.
!   - Wrap por eje solo si hay PBC en ese eje.
!   - Indexación de USS: idist = int(r*1000)+1 (clamp a [1,size(USS,1)]).
!   - DEVUELVE DELTV POSITIVO (energía de inserción adsorbato–adsorbato).
!─────────────────────────────────────────────────────────────────────
subroutine POTIN(MOLKIND, DELTV)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc, r2_min_image_frac
  use AdsorbateInput,  only : NMOLEC, NATOM, NATOMKIND
  use SimulationData,  only : RX, RY, RZ, RX1, RY1, RZ1, N, LOCATE, USS
  implicit none
  !---------------- argumentos
  integer,  intent(in)  :: MOLKIND
  real(rk), intent(out) :: DELTV
  !---------------- locales
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz
  integer  :: I, J, JIN, k1, k2, ipot, jpot
  real(rk) :: rix,riy,riz, rjx,rjy,rjz
  real(rk) :: s1,s2,s3, r2, rij
  integer  :: idist, maxr

  !── Celda "patrón" local: A/ACEL (misma escala que RX/RY/RZ reducidas)
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  DELTV = 0.0_rk
  maxr  = size(USS,1)   ! típicamente 5000

  !── Barrido por especies ya presentes
  do I = 1, NMOLEC
     if (N(I) == 0) cycle
     do J = 1, N(I)
        JIN = LOCATE(J, I)
        ! Interacción átomo-átomo entre la molécula de prueba (RX1/RY1/RZ1, MOLKIND)
        ! y la molécula existente (JIN,I)
        do k2 = 1, NATOM(I)
           jpot = NATOMKIND(k2, I)
           rjx  = RX(JIN, k2, I);  rjy = RY(JIN, k2, I);  rjz = RZ(JIN, k2, I)

           do k1 = 1, NATOM(MOLKIND)
              ipot = NATOMKIND(k1, MOLKIND)
              rix  = RX1(k1);       riy = RY1(k1);        riz = RZ1(k1)

              ! Δs = Ainv_R · (ri - rj), wrap según PBC
              s1 = Ainv(1,1)*(rix-rjx) + Ainv(1,2)*(riy-rjy) + Ainv(1,3)*(riz-rjz)
              s2 = Ainv(2,1)*(rix-rjx) + Ainv(2,2)*(riy-rjy) + Ainv(2,3)*(riz-rjz)
              s3 = Ainv(3,1)*(rix-rjx) + Ainv(3,2)*(riy-rjy) + Ainv(3,3)*(riz-rjz)
              call wrap_by_pbc(s1,s2,s3, px,py,pz)

              r2  = r2_min_image_frac(G, s1,s2,s3)
              rij = sqrt(r2)

              ! Lookup radial en USS
              idist = int(rij*1000.0_rk) + 1
              if (idist < 1)     idist = 1
              if (idist > maxr)  idist = maxr

              DELTV = DELTV + USS(idist, jpot, ipot)
           end do
        end do
     end do
  end do
end subroutine POTIN

