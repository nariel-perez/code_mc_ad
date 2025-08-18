!─────────────────────────────────────────────────────────────────────
! File: estadistica.f90
! Llena la matriz CNF(ICNF,JCNF,KCNF, especie, atomo) con contadores
! usando celda general y PBC parciales (triclinic-ready).
!─────────────────────────────────────────────────────────────────────
subroutine estadistica(ncellmat)
  use PBC_Mod,        only : rk, Cell_t => Cell
  use InputParams,    only : cell, ACEL
  use GeomUtils,      only : cell_to_metric, wrap_by_pbc
  use AdsorbateInput, only : NMOLEC, NATOM
  use SimulationData, only : RX,RY,RZ, LOCATE, N, CNF
  implicit none
  integer, intent(in) :: ncellmat

  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  logical  :: px,py,pz
  integer  :: i, j, k, ipull
  real(rk) :: x, y, z, s1, s2, s3
  integer  :: half, ic, jc, kc
  real(rk), parameter :: epsb = 1.0e-12_rk

  ! Celda patrón coherente con UADS/coords reducidas: A/ACEL
  cellR = cell
  cellR%A = cell%A / ACEL
  call cellR%update()
  call cell_to_metric(cellR, G, Ainv, px,py,pz)

  half = ncellmat/2   ! igual que tu legacy (rango [-half:half])

  do i = 1, NMOLEC
     if (N(i) <= 0) cycle
     do j = 1, N(i)
        ipull = LOCATE(j, i)
        do k = 1, NATOM(i)

           ! cartesiano reducido → fraccional
           x = RX(ipull,k,i);  y = RY(ipull,k,i);  z = RZ(ipull,k,i)
           s1 = Ainv(1,1)*x + Ainv(1,2)*y + Ainv(1,3)*z
           s2 = Ainv(2,1)*x + Ainv(2,2)*y + Ainv(2,3)*z
           s3 = Ainv(3,1)*x + Ainv(3,2)*y + Ainv(3,3)*z

           call wrap_by_pbc(s1,s2,s3, px,py,pz)

           ! si no hay PBC, exigir estar dentro del primario
           if (.not. px) then
              if ( (s1 < -0.5_rk) .or. (s1 >= 0.5_rk) ) cycle
           end if
           if (.not. py) then
              if ( (s2 < -0.5_rk) .or. (s2 >= 0.5_rk) ) cycle
           end if
           if (.not. pz) then
              if ( (s3 < -0.5_rk) .or. (s3 >= 0.5_rk) ) cycle
           end if

           ! Mapear fraccional → índice entero en [-half:half]
           ic = int( floor( (s1 + 0.5_rk) * real(ncellmat, rk) ) ) - half
           jc = int( floor( (s2 + 0.5_rk) * real(ncellmat, rk) ) ) - half
           kc = int( floor( (s3 + 0.5_rk) * real(ncellmat, rk) ) ) - half

           ! Clamp por si cae en el borde superior 0.5 - eps
           if (ic < -half) ic = -half; if (ic > half) ic = half
           if (jc < -half) jc = -half; if (jc > half) jc = half
           if (kc < -half) kc = -half; if (kc > half) kc = half

           CNF(ic, jc, kc, i, k) = CNF(ic, jc, kc, i, k) + 1.0_rk
        end do
     end do
  end do
end subroutine estadistica

