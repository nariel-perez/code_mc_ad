!────────────────────────────────────────────────────────────────────────
!  File: Potencial.f90   (versión celda genérica, Plan A)
!  Calcula UADS manteniendo la malla fija y los mismos nombres,
!  pero usando `min_image` sobre una celda ya *reducida* (Å → σ).
!────────────────────────────────────────────────────────────────────────
subroutine POTENCIAL(EPS, sigma, sigmetano, NC, RCUT, diel)
   !---------------- Módulos -------------------------------
   use InputParams,   only : cell, acel, acelx, acely, acelz, mat
   use PBC_Mod,       only : rk, Cell_t => Cell, min_image
   use PhysicalConstants, only : FCLEC, FACTORELEC
   use SimulationData
   implicit none

   !---------------- Argumentos ----------------------------
   real(rk), intent(in)    :: EPS, sigma, sigmetano, RCUT, diel
   integer,  intent(in)    :: NC

   !---------------- Variables locales ---------------------
   integer :: NKIND, INKIND, KINDI, IPOT, I, IJ, K, J
   real(rk) :: PI, RCELE
   real(rk) :: EPSIINKIND, SIGMINKIND, QINKIND
   real(rk) :: RXI, RYI, RZI, DELTV, SIGMA1, FACTOR, SIGSQ, SIGCUB
   real(rk) :: RMIN, RMINSQ, SR3, SR9, SR2, SR6, VLRC0, WLRC0, DELTW
   real(rk) :: RXIJ, RYIJ, RZIJ, RIJSQ, RIJ, VIJ, WIJ, VIJR, ZESACT, ZI
   real(rk) :: REDELEC
   integer  :: ios
   real(rk) :: xmax, ymax, zmax          ! aún impresos en resumen

   real(rk) :: dr(3)                     ! vector mínima-imagen

   !---------------- Celda reducida (Å → σ) ----------------
   type(Cell_t) :: cellR
   cellR = cell
   cellR%A = (cell%A / sigmetano) * sigma   ! escala
   call cellR%update()

   !---------------- Constantes ----------------------------
   PI    = acos(-1._rk)
   RCELE = 0.5_rk          ! corte electrostático (unidades σ)

   xmax = acelx / acel     ! sólo si quieres imprimirlos
   ymax = acely / acel
   zmax = acelz / acel

   !---------------- Leer LJ.dat ---------------------------
   open(unit=11, file='LJ.dat', status='old', action='read', iostat=ios)
   if (ios /= 0) stop 'ERROR: no se pudo abrir LJ.dat'
   read(11,*,iostat=ios) NKIND
   if (ios /= 0) stop 'ERROR leyendo NKIND'

   do INKIND = 1, NKIND
      read(11,*,iostat=ios) KINDI, EPSIINKIND, SIGMINKIND, QINKIND
      if (ios /= 0) stop 'ERROR leyendo parámetro LJ'
      EPSI(INKIND) = EPSIINKIND
      SIGM(INKIND) = SIGMINKIND
      Q(INKIND)    = QINKIND
      Q(INKIND)    = QINKIND * FACTORELEC * FCLEC
   end do
   close(11)

   !---------------- Tabla de potencial --------------------
   do IPOT = 1, NKIND
      write(*,*) 'POTENCIAL PARA ', IPOT,' de ', NKIND
      do I = -mat/2, mat/2
         RZI = real(I,rk) / real(mat,rk)
         do IJ = -mat/2, mat/2
            RYI = real(IJ,rk) / real(mat,rk)
            do K = -mat/2, mat/2
               RXI = real(K,rk) / real(mat,rk)
               VIJ = 0.0_rk
               VIJR = 0.0_rk
               DELTV = 0.0_rk
               DELTW = 0.0_rk

               do J = 1, NC
                  SIGMA1 = sigma * (SGC(J) + SIGM(IPOT)) / (2.0_rk*sigmetano)
                  FACTOR = sqrt(EPSI(IPOT)*EPSAC(J)) / EPS
                  SIGSQ  = SIGMA1*SIGMA1
                  SIGCUB = SIGSQ*SIGMA1
                  RMIN   = 0.5_rk*SIGMA1
                  RMINSQ = RMIN*RMIN
                  SR3 = (SIGMA1/RCUT)**3.0_rk
                  SR9 = SR3**3.0_rk

                  ! --------- mínima-imagen genérica ----------
                  dr = min_image(cellR, [RXI,RYI,RZI], [RXC(J),RYC(J),RZC(J)])
                  RXIJ = dr(1)
                  RYIJ = dr(2)
                  RZIJ = dr(3)
                  RIJSQ = sum(dr*dr)
                  ! ------------------------------------------

                  if (RIJSQ < RMINSQ) then
                     DELTV = 1.0e6_rk
                     exit
                  else
                     SR2 = SIGSQ / RIJSQ
                     SR6 = SR2*SR2*SR2
                     VIJ = SR6*(SR6 - 1.0_rk)*FACTOR
                     if (RIJSQ > (8.0_rk*RMIN)**2) VIJ = 0.0_rk

                     WIJ = SR6*(SR6 - 0.5_rk)*FACTOR

                     if (RIJSQ < RCELE**2) then
                        RIJ   = sqrt(RIJSQ)
                        ZESACT= QAC(J)
                        ZI = Q(IPOT)
                        VIJR  = ZI*ZESACT*(1.0_rk/RIJ - 1.0_rk/RCELE + (1.0_rk/RCELE**2)*(RIJ-RCELE))
                     else
                        VIJR = 0.0_rk
                     end if
                     DELTV = DELTV + 4.0_rk*VIJ + VIJR
                     DELTW = DELTW + WIJ
                  end if
               end do

               UADS(K,IJ,I,IPOT) = DELTV
               if (DELTV > 0.0_rk) DELTV = 0.0_rk
               DELTV = DELTV / EPS * 8.31_rk
            end do
         end do
         WRITE(*,*) DELTV,I,IPOT
      end do
      !write(*,*) 'FIN POTENCIAL PARA ', IPOT
   end do

   write(*,*) 'Energía calculada'
end subroutine POTENCIAL

