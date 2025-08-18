!─────────────────────────────────────────────────────────────────────
! File: Adpotin.f90    (Fortran 2003+)
! Propósito:
!   Energía de interacción adsorbato–superficie para la molécula de
!   PRUEBA (RX1/RY1/RZ1) de especie MOLKIND, usando la tabla UADS.
!
! Notas:
!   - Se trabaja en fraccionales de la "celda patrón" (cellR = A/ACEL).
!   - s = Ainv_R · r_atom ; wrap por eje sólo si hay PBC en ese eje.
!   - Indexación de UADS: k,ij,i = nint(s*mat) con clamp a [-m2:m2].
!   - ipot = NATOMKIND(k1, MOLKIND) (tipo atómico del átomo k1).
!─────────────────────────────────────────────────────────────────────
subroutine ADPOTIN(MOLKIND, DELTVA)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL, mat
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc
  use AdsorbateInput,  only : NATOM, NATOMKIND
  use SimulationData,  only : RX1, RY1, RZ1, UADS
  implicit none
  !---------------- argumentos
  integer,  intent(in)  :: MOLKIND
  real(rk), intent(out) :: DELTVA
  !---------------- locales
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz
  integer  :: k1, ipot
  real(rk) :: s1,s2,s3
  integer  :: kk, ijj, ii, m2

  ! Celda "patrón" local coherente con UADS (A/ACEL)
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  DELTVA = 0.0_rk
  m2     = mat/2

  ! Recorre átomos de la molécula de prueba
  do k1 = 1, NATOM(MOLKIND)
     ipot = NATOMKIND(k1, MOLKIND)

     ! s = Ainv_R · r_atom_trial  (r en unidades reducidas)
     s1 = Ainv(1,1)*RX1(k1) + Ainv(1,2)*RY1(k1) + Ainv(1,3)*RZ1(k1)
     s2 = Ainv(2,1)*RX1(k1) + Ainv(2,2)*RY1(k1) + Ainv(2,3)*RZ1(k1)
     s3 = Ainv(3,1)*RX1(k1) + Ainv(3,2)*RY1(k1) + Ainv(3,3)*RZ1(k1)

     ! Wrap por eje sólo si hay PBC en ese eje (e.g., BCZ=0 → no envolver z)
     call wrap_by_pbc(s1,s2,s3, px,py,pz)

     ! Índices de malla [-m2:m2] (K, IJ, I) con clamp en bordes
     kk  = nint( s1*real(mat,rk) );  if (kk  >  m2) kk  =  m2;  if (kk  < -m2) kk  = -m2
     ijj = nint( s2*real(mat,rk) );  if (ijj >  m2) ijj =  m2;  if (ijj < -m2) ijj = -m2
     ii  = nint( s3*real(mat,rk) );  if (ii  >  m2) ii  =  m2;  if (ii  < -m2) ii  = -m2

     ! Acumula potencial adsorbato–superficie desde la tabla UADS
     DELTVA = DELTVA + UADS(kk, ijj, ii, ipot)
  end do
end subroutine ADPOTIN

