!─────────────────────────────────────────────────────────────────────
! File: In.f90  (Fortran 2003+)
! Inserción GCMC con celda general y PBC parciales.
! Opción original: especie elegida UNIFORME y β = 1/T* (T reducido).
! Con chequeos defensivos para evitar segfaults.
!─────────────────────────────────────────────────────────────────────
subroutine IN(temp, Z, sigma, eps, rcut, V, VA, VG, W, CREATE, CR, jpasos, canonicalmolecules)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL, mat
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc
  use RotationModule,  only : GetRotationMatrix
  use AdsorbateInput,  only : NMOLEC, NATOM, RX0, RY0, RZ0
  use SimulationData,  only : RX1, RY1, RZ1, N
  implicit none
  !---------------- argumentos
  real(rk), intent(in)    :: temp, sigma, eps, rcut
  real(rk), intent(in)    :: Z(:)                 ! actividades por especie (ya precomputadas)
  real(rk), intent(inout) :: V, VA, VG, W
  logical,  intent(out)   :: CREATE
  real(rk), intent(out)   :: CR
  integer,  intent(in)    :: jpasos, canonicalmolecules
  !---------------- locales
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz, inside
  integer  :: molkind, k, ntrial, nAtomsSpec
  real(rk) :: u, beta, deltv, deltva, dU
  real(rk) :: s1,s2,s3, scom(3), rcom(3), R(3,3), dx,dy,dz
  real(rk), parameter :: pi   = 3.14159265358979323846_rk
  real(rk), parameter :: epsb = 1.0e-12_rk
  real(rk) :: DELTCB

  interface
     subroutine POTIN(molkind, deltv)
       use PBC_Mod, only: rk
       integer,  intent(in)  :: molkind
       real(rk), intent(out) :: deltv
     end subroutine
     subroutine ADPOTIN(molkind, deltv)
       use PBC_Mod, only: rk
       integer,  intent(in)  :: molkind
       real(rk), intent(out) :: deltv
     end subroutine
     subroutine ADD(molkind)
       integer, intent(in) :: molkind
     end subroutine
  end interface

  !---------------- seguridad mínima ----------------
  CREATE = .false.;  CR = 0.0_rk
  if (.not. allocated(NATOM)) return
  if (.not. allocated(RX1))   return
  if (.not. allocated(RY1))   return
  if (.not. allocated(RZ1))   return
  if (size(Z) < max(1,NMOLEC)) return

  !── Celda patrón coherente con tablas (A/ACEL)
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  !── β en UNIDADES REDUCIDAS (legacy): T* = T/ε  →  β = 1/T*
  if (temp <= 0.0_rk) return
  beta = 1.0_rk / temp

  !── Especie UNIFORME (opción original)
  call random_number(u)
  molkind = int( u * real(NMOLEC, rk) ) + 1
  if (molkind < 1) molkind = 1
  if (molkind > NMOLEC) molkind = NMOLEC

  nAtomsSpec = NATOM(molkind)
  if (nAtomsSpec <= 0) return
  ntrial = N(molkind) + 1

  !── Centro de masa trial en FRACCIONALES: U(-0.5,0.5) con wrap selectivo
  call random_number(scom(1))
  call random_number(scom(2))
  call random_number(scom(3))
  scom = scom - 0.5_rk
  call wrap_by_pbc(scom(1), scom(2), scom(3), px,py,pz)  ! no envuelve ejes sin PBC

  !── Orientación aleatoria (isotrópica aproximada)
  call random_number(u); dx = (u*2.0_rk - 1.0_rk)*pi
  call random_number(u); dy = (u*2.0_rk - 1.0_rk)*pi
  call random_number(u); dz = (u*2.0_rk - 1.0_rk)*pi
  call GetRotationMatrix(dx,dy,dz,R)

  !── Construir RX1/RY1/RZ1 = R*r0 + rCOM_trial (todo en reducidas)
  rcom = matmul(cellR%A, scom)
  do k = 1, nAtomsSpec
     RX1(k) = R(1,1)*RX0(k,molkind) + R(1,2)*RY0(k,molkind) + R(1,3)*RZ0(k,molkind) + rcom(1)
     RY1(k) = R(2,1)*RX0(k,molkind) + R(2,2)*RY0(k,molkind) + R(2,3)*RZ0(k,molkind) + rcom(2)
     RZ1(k) = R(3,1)*RX0(k,molkind) + R(3,2)*RY0(k,molkind) + R(3,3)*RZ0(k,molkind) + rcom(3)
  end do

  !── Si un eje NO tiene PBC, exigir que todos los átomos queden dentro del primario
  inside = .true.
  do k = 1, nAtomsSpec
     s1 = Ainv(1,1)*RX1(k) + Ainv(1,2)*RY1(k) + Ainv(1,3)*RZ1(k)
     s2 = Ainv(2,1)*RX1(k) + Ainv(2,2)*RY1(k) + Ainv(2,3)*RZ1(k)
     s3 = Ainv(3,1)*RX1(k) + Ainv(3,2)*RY1(k) + Ainv(3,3)*RZ1(k)
     if (.not. px) inside = inside .and. (abs(s1) <= 0.5_rk - epsb)
     if (.not. py) inside = inside .and. (abs(s2) <= 0.5_rk - epsb)
     if (.not. pz) inside = inside .and. (abs(s3) <= 0.5_rk - epsb)
     if (.not. inside) exit
  end do
  if (.not. inside) return

  !── Energía de inserción (adsorbato–superficie + adsorbato–adsorbato)
  !   Nota: si ADPOTIN/POTIN intentan indexar fuera de UADS, compilar con -fcheck=bounds.
  call ADPOTIN(molkind, deltva)
  if (.not. (deltva < 1.0e9_rk)) return        ! defensa por si vuelve NaN/Inf
  if (deltva > 1.0e3_rk) return                ! mismo criterio legacy

  call POTIN (molkind, deltv)
  if (.not. (deltv < 1.0e9_rk)) return

  dU = deltva + deltv

  !── Aceptación (legacy con DELTCB y corte 75)
  !   DELTCB = β ΔU - ln( Z / (N+1) )
  if (Z(molkind) <= 0.0_rk) return
  DELTCB = beta * dU - log( max(Z(molkind), 1.0e-300_rk) / real(ntrial, rk) )
  if (DELTCB < 75.0_rk) then
     call random_number(u)
     if ( (DELTCB <= 0.0_rk) .or. (exp(-DELTCB) > u) ) then
        call ADD(molkind)
        VG = VG + deltv
        VA = VA + deltva
        V  = VG + VA
        CREATE = .true.
        CR = 1.0_rk
     end if
  end if
end subroutine IN

