!─────────────────────────────────────────────────────────────────────
! File: Adpotout.f90   (Fortran 2003+)
! Propósito:
!   Energía adsorbato–superficie para la molécula EXISTENTE que será
!   removida (índice IPULL de la especie MOLKIND), usando la tabla UADS.
!
! Notas:
!   - Trabajo en fraccionales de la "celda patrón" (cellR = A/ACEL).
!   - s = Ainv_R · r_atom ; wrap por eje solo si hay PBC en ese eje.
!   - Indexación UADS: K, IJ, I = nint(s*mat) con clamp a [-m2:m2].
!   - Convención legacy de "out": DELTV = -DELTV al final.
!─────────────────────────────────────────────────────────────────────
subroutine ADPOTOUT(IPULL, MOLKIND, DELTV)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL, mat
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc
  use AdsorbateInput,  only : NATOM, NATOMKIND
  use SimulationData,  only : RX, RY, RZ, UADS
  implicit none
  !---------------- argumentos
  integer,  intent(in)  :: IPULL, MOLKIND
  real(rk), intent(out) :: DELTV
  !---------------- locales
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz
  integer  :: k1, ipot
  real(rk) :: s1,s2,s3
  integer  :: kk, ijj, ii, m2

  ! Celda "patrón" coherente con la construcción de UADS (A/ACEL)
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  DELTV = 0.0_rk
  m2    = mat/2

  ! Recorre átomos de la molécula a remover (coordenadas reducidas)
  do k1 = 1, NATOM(MOLKIND)
     ipot = NATOMKIND(k1, MOLKIND)

     ! s = Ainv_R · r_atom_existente
     s1 = Ainv(1,1)*RX(IPULL,k1,MOLKIND) + Ainv(1,2)*RY(IPULL,k1,MOLKIND) + Ainv(1,3)*RZ(IPULL,k1,MOLKIND)
     s2 = Ainv(2,1)*RX(IPULL,k1,MOLKIND) + Ainv(2,2)*RY(IPULL,k1,MOLKIND) + Ainv(2,3)*RZ(IPULL,k1,MOLKIND)
     s3 = Ainv(3,1)*RX(IPULL,k1,MOLKIND) + Ainv(3,2)*RY(IPULL,k1,MOLKIND) + Ainv(3,3)*RZ(IPULL,k1,MOLKIND)

     ! Wrap por eje solo donde hay PBC (ej.: BCZ=0 → NO envolver z)
     call wrap_by_pbc(s1,s2,s3, px,py,pz)

     ! Índices de malla [-m2:m2] (K, IJ, I) con clamp en bordes
     kk  = nint( s1*real(mat,rk) );  if (kk  >  m2) kk  =  m2;  if (kk  < -m2) kk  = -m2
     ijj = nint( s2*real(mat,rk) );  if (ijj >  m2) ijj =  m2;  if (ijj < -m2) ijj = -m2
     ii  = nint( s3*real(mat,rk) );  if (ii  >  m2) ii  =  m2;  if (ii  < -m2) ii  = -m2

     ! Acumula potencial adsorbato–superficie desde la tabla UADS
     DELTV = DELTV + UADS(kk, ijj, ii, ipot)
  end do

  ! Convención legacy para eliminación ("out")
  DELTV = -DELTV
end subroutine ADPOTOUT

