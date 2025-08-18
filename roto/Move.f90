!─────────────────────────────────────────────────────────────────
! File: Move.f90  (Fortran 2003 / rk kind)
! Movimiento MC (rotación/traslación) con:
!  - Selección uniforme de especie y molécula (estilo legacy)
!  - PBC genéricas (triclinic) vía A/Ainv y wrap selectivo por eje
!  - Salvavidas mínimos (sistema vacío, índices válidos)
!  - Convención de energías legacy:
!      POTOUT(ipull,kind,del)   → del = -E_old  (ff)
!      ADPOTOUT(ipull,kind,del) → del = -E_old  (fs)
!    y análogos para el trial
!─────────────────────────────────────────────────────────────────
subroutine MOVE(TEMP, Z, SIGMA, EPS, RCUT, V, VA, VG, W, GHOST, JPASOS)
  use PBC_Mod,         only : rk, Cell_t => Cell
  use InputParams,     only : cell, ACEL
  use GeomUtils,       only : cell_to_metric, wrap_by_pbc
  use RotationModule,  only : GetRotationMatrix
  use AdsorbateInput,  only : NMOLEC, NATOM, RX0, RY0, RZ0
  use SimulationData,  only : RX, RY, RZ, RX1, RY1, RZ1, N, LOCATE, &
                               ANX, ANGY, ANZ, EXNEW, EYNEW, EZNEW
  implicit none
  !---------------- args
  real(rk), intent(in)    :: TEMP, Z(:), SIGMA, EPS, RCUT
  real(rk), intent(inout) :: V, VA, VG, W
  logical,  intent(out)   :: GHOST
  integer,  intent(in)    :: JPASOS
  !---------------- locales
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz, inside
  integer  :: molkind, nloc, ipull, i, ij, nat
  real(rk) :: u, beta
  real(rk) :: VOLD, VGOLD, VAOLD, VANT, VGANT, VAANT
  real(rk) :: VNUEVA, VGNUEVA, VANUEVA
  real(rk) :: DELTVOUT, DELTVADOUT, DELTVOUT2, DELTVADOUT2
  real(rk) :: DELTCB
  real(rk) :: dx,dy,dz, R(3,3)
  real(rk) :: s1,s2,s3, epsb
  real(rk) :: deltax, deltays, deltaz, delta
  real(rk) :: rx_anchor, ry_anchor, rz_anchor
  real(rk) :: EXOLD, EYOLD, EZOLD
  real(rk), allocatable :: RXB(:), RYB(:), RZB(:)
  integer  :: totalN, sp

  real(rk), parameter :: pi   = 3.14159265358979323846_rk
  real(rk), parameter :: pi2  = 3.14159265358979323846_rk
  real(rk), parameter :: big  = 75.0_rk

  ! Las rutinas de energía se enlazan por nombre (estilo legacy)
  external :: POTOUT, ADPOTOUT

  !──────────────────────────────────────────────────────────────
  ! Salvavidas: si no hay partículas, nada que mover
  !──────────────────────────────────────────────────────────────
  totalN = 0
  do sp = 1, NMOLEC
     totalN = totalN + N(sp)
  end do
  if (totalN == 0) then
     GHOST = .false.
     return
  end if

  GHOST = .false.
  epsb  = 1.0e-12_rk

  !── Celda patrón coherente con UADS: A/ACEL
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  !── β reducido (legacy): β = 1/T*
  beta = 1.0_rk / max(TEMP, 1.0e-12_rk)

  !── Elegir especie UNIFORME y una molécula al azar de esa especie
  call random_number(u)
  molkind = int(u*real(NMOLEC,rk)) + 1
  if (molkind < 1) molkind = 1
  if (molkind > NMOLEC) molkind = NMOLEC
  if (N(molkind) <= 0) return

  call random_number(u)
  nloc  = int(u*real(N(molkind),rk)) + 1
  if (nloc < 1) nloc = 1
  if (nloc > N(molkind)) nloc = N(molkind)

  ipull = LOCATE(nloc, molkind)
  if (ipull < 1 .or. ipull > size(RX,1)) return

  nat   = NATOM(molkind)

  !── Copia de respaldo (viejas coords) en arrays locales
  allocate(RXB(nat), RYB(nat), RZB(nat))
  do i = 1, nat
     RXB(i) = RX(ipull,i,molkind)
     RYB(i) = RY(ipull,i,molkind)
     RZB(i) = RZ(ipull,i,molkind)
  end do
  EXOLD = ANX(ipull,molkind);  EYOLD = ANGY(ipull,molkind);  EZOLD = ANZ(ipull,molkind)

  !── Energías “salientes” (legacy: POTOUT/ADPOTOUT devuelven -E_old)
  VOLD = V; VGOLD = VG; VAOLD = VA
  call POTOUT  (ipull, molkind, DELTVOUT)
  call ADPOTOUT(ipull, molkind, DELTVADOUT)
  VANT  = VOLD  + DELTVOUT  + DELTVADOUT   ! = V - E_old + E_old?  (convención legacy)
  VGANT = VGOLD + DELTVOUT
  VAANT = VAOLD + DELTVADOUT

  !── Elegir ROTACIÓN (1) o TRASLACIÓN (2) — monatom: solo traslación
  call random_number(u)
  ij = int(u*2.0_rk) + 1
  if (nat == 1) ij = 2

  select case (ij)

  case (1)   !==================== ROTACIÓN ==========================
     ! Ancla legacy: átomo #2 si existe; si no, #1
     if (nat >= 2) then
        rx_anchor = RXB(2);  ry_anchor = RYB(2);  rz_anchor = RZB(2)
     else
        rx_anchor = RXB(1);  ry_anchor = RYB(1);  rz_anchor = RZB(1)
     end if

     ! Ángulos uniformes en [-π, π]
     call random_number(u); dx = (2.0_rk*u - 1.0_rk)*pi2
     call random_number(u); dy = (2.0_rk*u - 1.0_rk)*pi2
     call random_number(u); dz = (2.0_rk*u - 1.0_rk)*pi2
     call GetRotationMatrix(dx,dy,dz,R)

     ! r_trial = R*r0 + r_anchor  (usamos molde RX0/RY0/RZ0)
     do i = 1, nat
        RX1(i) = R(1,1)*RX0(i,molkind) + R(1,2)*RY0(i,molkind) + R(1,3)*RZ0(i,molkind) + rx_anchor
        RY1(i) = R(2,1)*RX0(i,molkind) + R(2,2)*RY0(i,molkind) + R(2,3)*RZ0(i,molkind) + ry_anchor
        RZ1(i) = R(3,1)*RX0(i,molkind) + R(3,2)*RY0(i,molkind) + R(3,3)*RZ0(i,molkind) + rz_anchor
     end do

     ! Wrap selectivo + “inside” si un eje no tiene PBC; volver a cart trial
     inside = .true.
     do i = 1, nat
        s1 = Ainv(1,1)*RX1(i) + Ainv(1,2)*RY1(i) + Ainv(1,3)*RZ1(i)
        s2 = Ainv(2,1)*RX1(i) + Ainv(2,2)*RY1(i) + Ainv(2,3)*RZ1(i)
        s3 = Ainv(3,1)*RX1(i) + Ainv(3,2)*RY1(i) + Ainv(3,3)*RZ1(i)
        call wrap_by_pbc(s1,s2,s3, px,py,pz)
        if (.not. px) inside = inside .and. (abs(s1) <= 0.5_rk - epsb)
        if (.not. py) inside = inside .and. (abs(s2) <= 0.5_rk - epsb)
        if (.not. pz) inside = inside .and. (abs(s3) <= 0.5_rk - epsb)
        if (.not. inside) then
           call restore_and_exit()
           return
        end if
        RX1(i) = cellR%A(1,1)*s1 + cellR%A(1,2)*s2 + cellR%A(1,3)*s3
        RY1(i) = cellR%A(2,1)*s1 + cellR%A(2,2)*s2 + cellR%A(2,3)*s3
        RZ1(i) = cellR%A(3,1)*s1 + cellR%A(3,2)*s2 + cellR%A(3,3)*s3
     end do
     EXNEW = dx/pi2; EYNEW = dy/pi2; EZNEW = dz/pi2

  case (2)   !==================== TRASLACIÓN ========================
     delta  = 0.01_rk
     call random_number(u); deltax = (u - 0.5_rk) * SIGMA * delta
     call random_number(u); deltays= (u - 0.5_rk) * SIGMA * delta
     call random_number(u); deltaz = (u - 0.5_rk) * SIGMA * delta

     do i = 1, nat
        RX1(i) = RXB(i) + deltax
        RY1(i) = RYB(i) + deltays
        RZ1(i) = RZB(i) + deltaz
     end do

     inside = .true.
     do i = 1, nat
        s1 = Ainv(1,1)*RX1(i) + Ainv(1,2)*RY1(i) + Ainv(1,3)*RZ1(i)
        s2 = Ainv(2,1)*RX1(i) + Ainv(2,2)*RY1(i) + Ainv(2,3)*RZ1(i)
        s3 = Ainv(3,1)*RX1(i) + Ainv(3,2)*RY1(i) + Ainv(3,3)*RZ1(i)
        call wrap_by_pbc(s1,s2,s3, px,py,pz)
        if (.not. px) inside = inside .and. (abs(s1) <= 0.5_rk - epsb)
        if (.not. py) inside = inside .and. (abs(s2) <= 0.5_rk - epsb)
        if (.not. pz) inside = inside .and. (abs(s3) <= 0.5_rk - epsb)
        if (.not. inside) then
           call restore_and_exit()
           return
        end if
        RX1(i) = cellR%A(1,1)*s1 + cellR%A(1,2)*s2 + cellR%A(1,3)*s3
        RY1(i) = cellR%A(2,1)*s1 + cellR%A(2,2)*s2 + cellR%A(2,3)*s3
        RZ1(i) = cellR%A(3,1)*s1 + cellR%A(3,2)*s2 + cellR%A(3,3)*s3
     end do
     EXNEW = EXOLD; EYNEW = EYOLD; EZNEW = EZOLD
  end select

  !── Escribir TRIAL en el estado y evaluar energías nuevas
  do i = 1, nat
     RX(ipull,i,molkind) = RX1(i)
     RY(ipull,i,molkind) = RY1(i)
     RZ(ipull,i,molkind) = RZ1(i)
  end do
  ANX(ipull,molkind) = EXNEW;  ANGY(ipull,molkind) = EYNEW;  ANZ(ipull,molkind) = EZNEW

  call POTOUT  (ipull, molkind, DELTVOUT2)     ! = -E_trial(ff)
  call ADPOTOUT(ipull, molkind, DELTVADOUT2)   ! = -E_trial(fs)

  VNUEVA  = VANT  - DELTVOUT2  - DELTVADOUT2   ! = (V - E_old) + E_trial
  VGNUEVA = VGANT - DELTVOUT2
  VANUEVA = VAANT - DELTVADOUT2

  DELTCB = beta * (VNUEVA - VOLD)

  if (DELTCB < big) then
     call random_number(u)
     if ( (DELTCB <= 0.0_rk) .or. (exp(-DELTCB) > u) ) then
        V  = VNUEVA
        VG = VGNUEVA
        VA = VANUEVA
        if (allocated(RXB)) deallocate(RXB,RYB,RZB)
        return
     else
        call restore_and_exit()
        return
     end if
  else
     call restore_and_exit()
     return
  end if

contains

  subroutine restore_and_exit()
    integer :: j
    do j = 1, nat
       RX(ipull,j,molkind) = RXB(j)
       RY(ipull,j,molkind) = RYB(j)
       RZ(ipull,j,molkind) = RZB(j)
    end do
    ANX(ipull,molkind) = EXOLD
    ANGY(ipull,molkind) = EYOLD
    ANZ(ipull,molkind) = EZOLD
    V  = VOLD; VG = VGOLD; VA = VAOLD
    if (allocated(RXB)) deallocate(RXB,RYB,RZB)
  end subroutine restore_and_exit

end subroutine MOVE

