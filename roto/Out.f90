!─────────────────────────────────────────────────────────────────
! File: Out.f90   (Fortran 2003)
! Lógica clásica: elegir especie, elegir molécula, calcular ΔU,
! aceptar con: arg = β*ΔU − log( N(isp) / Z(isp) )
! Salvavidas mínimos:
!   - si el sistema está vacío → ghost=.true.; return
!   - si la especie elegida no tiene población → ghost=.true.; return
!   - si el índice/slot son inválidos → ghost=.true.; return
! No se definen INTERFACEs aquí (implícito como en el legacy).
!─────────────────────────────────────────────────────────────────
subroutine OUT(tempstar, Z, sigma, eps, rcut, V, VA, VG, W, ghost, jpas, NMIN, NMAXI, canonicalmolecules)
   use PBC_Mod,         only : rk
   use SimulationData,  only : RX, RY, RZ, UADS, USS, LOCATE, N
   use AdsorbateInput,  only : NMOLEC, NATOM
   implicit none
   !--------------- argumentos -----------------
   real(rk), intent(in)    :: tempstar, sigma, eps, rcut
   real(rk), intent(inout) :: V, VA, VG, W
   real(rk), intent(in)    :: Z(:)
   logical, intent(out)    :: ghost
   integer,  intent(in)    :: jpas
   integer,  intent(in)    :: NMIN(:), NMAXI(:)
   integer,  intent(in)    :: canonicalmolecules
   !--------------- locales --------------------
   integer  :: i, isp, j, ipull, ntot
   real(rk) :: beta, dU, deltv, deltva
   real(rk) :: u, ur, zsafe, arg

   ! (En el legacy, POTOUT/ADPOTOUT/REMOVE se enlazan por nombre;
   !  no hace falta interface aquí si las firmas coinciden.)
   external :: POTOUT, ADPOTOUT, REMOVE

   ghost = .false.
   beta  = 1.0_rk / max(tempstar, 1.0e-12_rk)

   ! ===== salvavidas mínimo #0: ¿hay algo para remover? =====
   ntot = 0
   do i = 1, NMOLEC
      ntot = ntot + N(i)
   end do
   if (ntot == 0) then
      ghost = .true.
      return
   end if

   ! ===== elegir especie uniformemente (estilo legacy) =====
   call random_number(ur)
   isp = int(ur*real(size(N),rk)) + 1
   if (isp < 1) isp = 1
   if (isp > size(N)) isp = size(N)

   ! si no hay población en esa especie, abortar silenciosamente (legacy-like)
   if (N(isp) <= 0) then
      ghost = .true.
      return
   end if

   ! (si usás límites canónicos por especie, respetalos)
   if (N(isp) <= NMIN(isp)) then
      ghost = .true.
      return
   end if

   ! ===== elegir la j-ésima molécula de esa especie =====
   call random_number(ur)
   j = int(ur*real(N(isp), rk)) + 1
   if (j < 1 .or. j > N(isp)) then
      ghost = .true.
      return
   end if

   ! mapear al "slot" físico
   ipull = LOCATE(j, isp)
   if (ipull < 1 .or. ipull > size(RX,1)) then
      ghost = .true.
      return
   end if

   ! ===== energía del removido (ΔU = ΔV_ads-surf + ΔV_ads-ads) =====
   call POTOUT  (ipull, isp, deltv)     ! parte adsorbato-superficie
   call ADPOTOUT(ipull, isp, deltva)    ! parte adsorbato-adsorbato
   dU = deltv + deltva                  ! (en unidades consistentes con el loop, típicamente ε)

   ! ===== aceptación de borrado (formulación estándar legacy) =====
   zsafe = max(Z(isp), 1.0e-300_rk)
   arg   = beta*dU - log( real(max(N(isp),1), rk) / zsafe )
   call random_number(u)

   if (arg <= 0.0_rk .or. exp(-arg) > u) then
      ! aceptado → actualizar estado y energías (en las mismas unidades del loop)
      call REMOVE(isp, ipull)
      V  = V  - dU
      VA = VA - deltva
      VG = VG - (dU - deltva)
   else
      ghost = .true.
   end if
end subroutine OUT

