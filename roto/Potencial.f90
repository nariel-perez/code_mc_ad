!------------------------------------------------------------------------
! POTENCIAL (UADS) — versión rápida, general (triclinic, sin OpenMP)
! - Fraccionales + métrica G (celda patrón local cellR)
! - Wrap por eje según PBC
! - Precálculos por IPOT
! - Evita matmul/sqrt en el núcleo
!------------------------------------------------------------------------
subroutine POTENCIAL(EPS, sigma, sigmetano, NC, RCUT, diel)
  use InputParams,        only: mat, cell                ! 'cell' con a,b,c,angulos, pbc
  use PBC_Mod,            only: rk, Cell_t => Cell       ! tipo y kind
  use PhysicalConstants,  only: FCLEC, FACTORELEC
  use SimulationData,     only: UADS, RXC, RYC, RZC, EPSAC, SGC, EPSI, SIGM, Q, QAC
  implicit none

  !------------------------- argumentos
  integer, intent(in) :: NC
  real(rk), intent(in) :: EPS, sigma, sigmetano, RCUT, diel

  !------------------------- locales
  integer :: NKIND, INKIND, KINDI, IPOT
  integer :: i, ij, k, j, ios, m2
  real(rk) :: sI1, sI2, sI3, s1, s2, s3, r2, rij
  real(rk) :: sr2, sr6, vij, wij, vijr, deltv, deltw, last_deltv,pp
  logical :: px, py, pz
  real(rk), parameter :: RCELE = 0.5_rk
  real(rk), parameter :: BIGV  = 1.0e6_rk
  real(rk), parameter :: ZERO  = 0.0_rk

  ! Celda patrón local (escala legacy) y métricas
  type(Cell_t) :: cellR
  real(rk) :: G(3,3), Ainv(3,3)
  real(rk) :: g11,g22,g33,g12,g13,g23

  ! Sólido en fraccionales (celda patrón)
  real(rk), allocatable :: sC1(:), sC2(:), sC3(:)

  ! Precálculos por átomo (dependen de IPOT)
  real(rk), allocatable :: sigma1(:), sigsq(:), rmin(:), r2min(:), r2cut(:), factor(:)

  !------------------------- LJ.dat (legacy: primera línea = NKIND)
  open(11, file='LJ.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,*) 'Error al abrir LJ.dat'
     return
  end if
  read(11,*,iostat=ios) NKIND
  if (ios /= 0 .or. NKIND <= 0) then
     write(*,*) 'Error al leer NKIND en LJ.dat'
     close(11); return
  end if
  do INKIND = 1, NKIND
     read(11,*,iostat=ios) KINDI, EPSI(INKIND), SIGM(INKIND), Q(INKIND)
     if (ios /= 0) then
        write(*,*) 'Error al leer línea ', INKIND, ' de LJ.dat'
        close(11); return
     end if
     Q(INKIND) = Q(INKIND) * FACTORELEC
     Q(INKIND) = Q(INKIND) * FCLEC
  end do
  close(11)

  !------------------------- Celda patrón local (Å -> reducidas legacy)
  cellR = cell
  ! equivalentes: (cell%A / sigmetano) * sigma == cell%A / ACEL   (ACEL = sigmetano/ sigma^-1)
  cellR%A = (cell%A / sigmetano) * sigma
  call cellR%update()   ! debe dejar Ainv listo

  Ainv = cellR%Ainv
  G    = transpose(cellR%A)  ! G = A^T * A
  G    = matmul(G, cellR%A)

  g11 = G(1,1); g22 = G(2,2); g33 = G(3,3)
  g12 = G(1,2); g13 = G(1,3); g23 = G(2,3)

  px = cellR%pbc(1); py = cellR%pbc(2); pz = cellR%pbc(3)

  !------------------------- Sólido en fraccionales (una vez)
  allocate(sC1(NC), sC2(NC), sC3(NC))
  do j = 1, NC
     s1 = Ainv(1,1)*RXC(j) + Ainv(1,2)*RYC(j) + Ainv(1,3)*RZC(j)
     s2 = Ainv(2,1)*RXC(j) + Ainv(2,2)*RYC(j) + Ainv(2,3)*RZC(j)
     s3 = Ainv(3,1)*RXC(j) + Ainv(3,2)*RYC(j) + Ainv(3,3)*RZC(j)
     sC1(j) = s1;  sC2(j) = s2;  sC3(j) = s3
  end do

  !------------------------- Malla
  m2 = mat/2
  write(*,*) '-----------------------'

  do IPOT = 1, NKIND
     write(*,*) 'POTENCIAL PARA ', IPOT, ' de ', NKIND

     !------ precálculos por átomo del sólido (dependen de IPOT)
     allocate(sigma1(NC), sigsq(NC), rmin(NC), r2min(NC), r2cut(NC), factor(NC))
     do j = 1, NC
        sigma1(j) = sigma * ( SGC(j) + SIGM(IPOT) ) / ( 2.0_rk*sigmetano )
        sigsq(j)  = sigma1(j)*sigma1(j)
        rmin(j)   = 0.5_rk*sigma1(j)
        r2min(j)  = rmin(j)*rmin(j)
        r2cut(j)  = (8.0_rk*rmin(j))*(8.0_rk*rmin(j))
        factor(j) = sqrt( EPSI(IPOT)*EPSAC(j) ) / EPS
     end do

     last_deltv = ZERO

     !------ barrido de malla (fraccional)
     do i  = -m2, m2
        sI3 = real(i, rk)/real(mat, rk)
        do ij = -m2, m2
           sI2 = real(ij, rk)/real(mat, rk)
           do k = -m2, m2
              sI1 = real(k, rk)/real(mat, rk)

              deltv = ZERO
              deltw = ZERO

              do j = 1, NC
                 ! Δs con wrap por eje según PBC
                 s1 = sI1 - sC1(j);  if (px) s1 = s1 - nint(s1)
                 s2 = sI2 - sC2(j);  if (py) s2 = s2 - nint(s2)
                 s3 = sI3 - sC3(j);  if (pz) s3 = s3 - nint(s3)

                 ! r^2 = s^T G s
                 r2 = g11*s1*s1 + g22*s2*s2 + g33*s3*s3 + &
                       2.0_rk*(g12*s1*s2 + g13*s1*s3 + g23*s2*s3)

                 if (r2 < r2min(j)) then
                    deltv = BIGV
                    exit
                 else
                    sr2 = sigsq(j)/r2
                    sr6 = sr2*sr2*sr2
                    vij = sr6*(sr6 - 1.0_rk)*factor(j)
                    if (r2 > r2cut(j)) vij = ZERO

                    wij = sr6*(sr6 - 0.5_rk)*factor(j)

                    vijr = ZERO
                    if (r2 < RCELE*RCELE) then
                       rij  = sqrt(r2)
                       vijr = Q(IPOT)*QAC(j) * ( 1.0_rk/rij - 1.0_rk/RCELE + (1.0_rk/(RCELE*RCELE))*(rij - RCELE) )
                    end if
                 end if

                 deltv = deltv + 4.0_rk*vij + vijr
                 deltw = deltw + wij
              end do

              UADS(k, ij, i, IPOT) = deltv
              pp = deltv
              if (pp > 0.0_rk) pp = 0.0_rk
              pp = pp / EPS * 8.31_rk
              if (k == m2 .and. ij == m2) last_deltv = pp
           end do
        end do
        !write(*,*) last_deltv, i, IPOT
     end do

     deallocate(sigma1, sigsq, rmin, r2min, r2cut, factor)
  end do

  write(*,*) 'Energia Calculada'
  return
end subroutine POTENCIAL

