!─────────────────────────────────────────────────────────────────────
! File: In.f90  (Fortran 2003+)
! Inserción GCMC con celda general y PBC parciales.
! Opción original: especie elegida UNIFORME y β = 1/T* (T reducido).
! Incluye:
!  - Chequeo “inside” en ejes sin PBC
!  - Micro-ajuste opcional del COM en ejes sin PBC (evita rechazos triviales)
!  - Contadores de diagnóstico (try/inside/core/acc)
!─────────────────────────────────────────────────────────────────────
subroutine IN(temp, Z, sigma, eps, rcut, V, VA, VG, W, CREATE, CR, jpasos, canonicalmolecules)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc
  use RotationModule,  only : GetRotationMatrix
  use AdsorbateInput,  only : NMOLEC, NATOM, NATOMKIND, RX0, RY0, RZ0
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
  integer  :: molkind, k, ipot, ntrial
  real(rk) :: u, beta, deltv, deltva, dU
  real(rk) :: s1,s2,s3, scom(3), rcom(3), R(3,3), dx,dy,dz
  real(rk), parameter :: pi   = 3.14159265358979323846_rk
  real(rk), parameter :: epsb = 1.0e-12_rk
  real(rk) :: DELTCB
  !— Diagnóstico (persisten entre llamadas)
  integer, save :: ins_try=0, ins_inside_fail=0, ins_core=0, ins_acc=0
  !— Toggle del micro-ajuste del COM (para legacy exacto, ponelo en .false.)
  logical,  parameter :: do_adjust_com = .true.

  integer :: t
  real(rk) :: smin, smax, shift

  
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

  CREATE = .false.;  CR = 0.0_rk

  !── Celda patrón coherente con las tablas (A/ACEL)
  cellR   = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  !── β en UNIDADES REDUCIDAS (legacy): T* = T/ε  →  β = 1/T*
  beta = 1.0_rk / temp

  !── Especie UNIFORME (opción original)
  call random_number(u)
  molkind = int( u * real(NMOLEC, rk) ) + 1
  if (molkind < 1) molkind = 1
  if (molkind > NMOLEC) molkind = NMOLEC
  ntrial = N(molkind) + 1

  !── Centro de masa trial en FRACCIONALES: U(-0.5,0.5) por eje
  call random_number(scom(1))
  call random_number(scom(2))
  call random_number(scom(3))
  scom = scom - 0.5_rk
  call wrap_by_pbc(scom(1), scom(2), scom(3), px,py,pz)  ! no envuelve ejes sin PBC

  !── Orientación aleatoria (isotrópica aproximada con Euler independiente)
  call random_number(u); dx = (u*2.0_rk - 1.0_rk)*pi
  call random_number(u); dy = (u*2.0_rk - 1.0_rk)*pi
  call random_number(u); dz = (u*2.0_rk - 1.0_rk)*pi
  call GetRotationMatrix(dx,dy,dz,R)

  !── Construir RX1/RY1/RZ1 = R*r0 + rCOM_trial (todo en reducidas, celda A/ACEL)
  rcom = matmul(cellR%A, scom)
  do k = 1, NATOM(molkind)
     ipot  = NATOMKIND(k, molkind)
     RX1(k) = R(1,1)*RX0(k,molkind) + R(1,2)*RY0(k,molkind) + R(1,3)*RZ0(k,molkind) + rcom(1)
     RY1(k) = R(2,1)*RX0(k,molkind) + R(2,2)*RY0(k,molkind) + R(2,3)*RZ0(k,molkind) + rcom(2)
     RZ1(k) = R(3,1)*RX0(k,molkind) + R(3,2)*RY0(k,molkind) + R(3,3)*RZ0(k,molkind) + rcom(3)
  end do

  ins_try = ins_try + 1

  !── Si un eje NO tiene PBC, exigir que todos los átomos queden dentro del primario
  inside = .true.
  do k = 1, NATOM(molkind)
     s1 = Ainv(1,1)*RX1(k) + Ainv(1,2)*RY1(k) + Ainv(1,3)*RZ1(k)
     s2 = Ainv(2,1)*RX1(k) + Ainv(2,2)*RY1(k) + Ainv(2,3)*RZ1(k)
     s3 = Ainv(3,1)*RX1(k) + Ainv(3,2)*RY1(k) + Ainv(3,3)*RZ1(k)
     if (.not. px) inside = inside .and. (abs(s1) <= 0.5_rk - epsb)
     if (.not. py) inside = inside .and. (abs(s2) <= 0.5_rk - epsb)
     if (.not. pz) inside = inside .and. (abs(s3) <= 0.5_rk - epsb)
     if (.not. inside) exit
  end do

  if (.not. inside) then
     if (do_adjust_com) then
        !— Micro-ajuste del COM por cada eje sin PBC para meter el cluster en la celda
        

        do t = 1, 3
           ! ¿Este eje carece de PBC?
           if ( (t==1 .and. .not. px) .or. (t==2 .and. .not. py) .or. (t==3 .and. .not. pz) ) then
              smin = +1.0e99_rk
              smax = -1.0e99_rk

              do k = 1, NATOM(molkind)
                 s1 = Ainv(1,1)*RX1(k) + Ainv(1,2)*RY1(k) + Ainv(1,3)*RZ1(k)
                 s2 = Ainv(2,1)*RX1(k) + Ainv(2,2)*RY1(k) + Ainv(2,3)*RZ1(k)
                 s3 = Ainv(3,1)*RX1(k) + Ainv(3,2)*RY1(k) + Ainv(3,3)*RZ1(k)
                 if (t==1) then
                    smin = min(smin, s1);  smax = max(smax, s1)
                 elseif (t==2) then
                    smin = min(smin, s2);  smax = max(smax, s2)
                 else
                    smin = min(smin, s3);  smax = max(smax, s3)
                 end if
              end do

              if (smin < -0.5_rk .or. smax > 0.5_rk) then
                 if (smin < -0.5_rk) shift = (-0.5_rk - smin)
                 if (smax >  0.5_rk) shift = ( 0.5_rk - smax)

                 do k = 1, NATOM(molkind)
                    select case (t)
                    case (1)
                       RX1(k) = RX1(k) + shift*cellR%A(1,1)
                       RY1(k) = RY1(k) + shift*cellR%A(2,1)
                       RZ1(k) = RZ1(k) + shift*cellR%A(3,1)
                    case (2)
                       RX1(k) = RX1(k) + shift*cellR%A(1,2)
                       RY1(k) = RY1(k) + shift*cellR%A(2,2)
                       RZ1(k) = RZ1(k) + shift*cellR%A(3,2)
                    case (3)
                       RX1(k) = RX1(k) + shift*cellR%A(1,3)
                       RY1(k) = RY1(k) + shift*cellR%A(2,3)
                       RZ1(k) = RZ1(k) + shift*cellR%A(3,3)
                    end select
                 end do
              end if
           end if
        end do

        !— Rechequear inside una vez
        inside = .true.
        do k = 1, NATOM(molkind)
           s1 = Ainv(1,1)*RX1(k) + Ainv(1,2)*RY1(k) + Ainv(1,3)*RZ1(k)
           s2 = Ainv(2,1)*RX1(k) + Ainv(2,2)*RY1(k) + Ainv(2,3)*RZ1(k)
           s3 = Ainv(3,1)*RX1(k) + Ainv(3,2)*RY1(k) + Ainv(3,3)*RZ1(k)
           if (.not. px) inside = inside .and. (abs(s1) <= 0.5_rk - epsb)
           if (.not. py) inside = inside .and. (abs(s2) <= 0.5_rk - epsb)
           if (.not. pz) inside = inside .and. (abs(s3) <= 0.5_rk - epsb)
           if (.not. inside) exit
        end do
     end if

     if (.not. inside) then
        ins_inside_fail = ins_inside_fail + 1
        return
     end if
  end if

  !── Energía de inserción (adsorbato–superficie + adsorbato–adsorbato)
  call ADPOTIN(molkind, deltva)
  if (deltva > 1.0e5_rk) then
     ins_core = ins_core + 1
     return   ! criterio legacy: core duro → descartar
  end if

  call POTIN (molkind, deltv)
  dU = deltva + deltv

  !── Aceptación (formulación legacy con DELTCB y corte 75)
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
        ins_acc = ins_acc + 1
     end if
  end if

  !— Estadísticas de inserción cada 500 pasos de jpasos (tuneable)
  if (mod(jpasos, 500) == 0) then
     write(*,'(A,4I10)') 'IN stats:  try, inside_rej, core_rej, acc = ', &
          ins_try, ins_inside_fail, ins_core, ins_acc
  end if
end subroutine IN

