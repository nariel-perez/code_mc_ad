!─────────────────────────────────────────────────────────────────
! File: Adpotin.f90  (triclinic + PBC parciales, opción original)
!─────────────────────────────────────────────────────────────────
subroutine ADPOTIN(molkind, deltv)
  use PBC_Mod,        only : rk, Cell_t => Cell
  use InputParams,    only : cell, ACEL, mat
  use GeomUtils,      only : cell_to_metric, wrap_by_pbc
  use AdsorbateInput, only : NATOM, NATOMKIND
  use SimulationData, only : RX1, RY1, RZ1, UADS
  implicit none
  integer,  intent(in)  :: molkind
  real(rk), intent(out) :: deltv

  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz
  integer  :: k, ipot, ix, iy, iz, aux
  real(rk) :: s1, s2, s3, e

  ! Celda patrón (idéntica a la usada en POTENCIAL)
  cellR   = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  aux   = mat/2
  deltv = 0.0_rk

  do k = 1, NATOM(molkind)
     ipot = NATOMKIND(k, molkind)

     ! frac = Ainv * r_trial
     s1 = Ainv(1,1)*RX1(k) + Ainv(1,2)*RY1(k) + Ainv(1,3)*RZ1(k)
     s2 = Ainv(2,1)*RX1(k) + Ainv(2,2)*RY1(k) + Ainv(2,3)*RZ1(k)
     s3 = Ainv(3,1)*RX1(k) + Ainv(3,2)*RY1(k) + Ainv(3,3)*RZ1(k)

     ! Envolver SOLO en ejes con PBC (coherente con tabla y 2-D)
     call wrap_by_pbc(s1, s2, s3, px,py,pz)

     ! Índices de malla coherentes con POTENCIAL: (K, IJ, I) ≡ (x,y,z)
     ix = nint( s1 * real(mat, rk) )
     iy = nint( s2 * real(mat, rk) )
     iz = nint( s3 * real(mat, rk) )

     ! Clamps de seguridad
     if (ix < -aux) ix = -aux
     if (ix >  aux) ix =  aux
     if (iy < -aux) iy = -aux
     if (iy >  aux) iy =  aux
     if (iz < -aux) iz = -aux
     if (iz >  aux) iz =  aux

     e = UADS(ix, iy, iz, ipot)

     ! Si cayó en el core (sentinela grande), rechaza inserción
     if (e > 1.0e5_rk) then
        deltv = 1.0e6_rk
        return
     end if

     deltv = deltv + e
  end do
end subroutine ADPOTIN

