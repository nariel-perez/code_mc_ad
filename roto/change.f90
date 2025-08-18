!─────────────────────────────────────────────────────────────────────
! File: change.f90
! Stub “legacy-compatible”: no cambia parámetros ni el estado.
! Deja el gancho por si luego querés ajustar amplitudes/adaptativos.
!─────────────────────────────────────────────────────────────────────
subroutine change(temp, Z, sigma, eps, rcut, V, VA, VG, W, GHOST, jpasos)
  use PBC_Mod, only : rk
  implicit none
  real(rk), intent(in)    :: temp, sigma, eps, rcut
  real(rk), intent(in)    :: Z(:)
  real(rk), intent(inout) :: V, VA, VG, W
  logical,  intent(out)   :: GHOST
  integer,  intent(in)    :: jpasos

  ! Por compatibilidad con el flujo legacy
  GHOST = .false.
  ! (No-op)
end subroutine change

