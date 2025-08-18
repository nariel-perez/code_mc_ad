module GeomUtils
  use PBC_Mod, only : rk, Cell_t => Cell
  implicit none
  private
  public :: cell_to_metric, r2_min_image_frac, wrap_by_pbc

contains

  pure subroutine cell_to_metric(c, G, Ainv, pbcx, pbcy, pbcz)
    type(Cell_t), intent(in)  :: c
    real(rk),     intent(out) :: G(3,3), Ainv(3,3)
    logical,      intent(out) :: pbcx, pbcy, pbcz
    Ainv = c%Ainv
    G    = matmul(transpose(c%A), c%A)
    pbcx = c%pbc(1);  pbcy = c%pbc(2);  pbcz = c%pbc(3)
  end subroutine cell_to_metric

  pure real(rk) function r2_min_image_frac(G, s1, s2, s3) result(r2)
    real(rk), intent(in) :: G(3,3), s1, s2, s3
    r2 = G(1,1)*s1*s1 + G(2,2)*s2*s2 + G(3,3)*s3*s3 &
       + 2.0_rk*( G(1,2)*s1*s2 + G(1,3)*s1*s3 + G(2,3)*s2*s3 )
  end function r2_min_image_frac

  pure subroutine wrap_by_pbc(s1,s2,s3,pbcx,pbcy,pbcz)
    real(rk), intent(inout) :: s1,s2,s3
    logical,  intent(in)    :: pbcx,pbcy,pbcz
    if (pbcx) s1 = s1 - nint(s1)
    if (pbcy) s2 = s2 - nint(s2)
    if (pbcz) s3 = s3 - nint(s3)
  end subroutine wrap_by_pbc

end module GeomUtils

