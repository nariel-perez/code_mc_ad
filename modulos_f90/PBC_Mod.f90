!─────────────────────────────────────────────────────────────────────
!  File: PBC_Mod.f90                       (Fortran 2003 / 2008)
!  Manejo genérico de celdas periódicas: construcción, conversiones
!  cart ↔ frac, envolturas y vector mínima-imagen.
!─────────────────────────────────────────────────────────────────────
module PBC_Mod
   implicit none
   private                    ! todo oculto; se exporta al final

   !============================= KIND ===============================
   integer, parameter,public :: rk = selected_real_kind(12, 99)

   !============================= TYPE ===============================
   type, public :: Cell
      ! --- datos geométricos ---------------------------------------
      real(rk) :: a_len = 0._rk, b_len = 0._rk, c_len = 0._rk
      real(rk) :: alpha_deg = 90._rk, beta_deg = 90._rk, gamma_deg = 90._rk
      ! vectores de red en columnas: A(:,1)=a, A(:,2)=b, A(:,3)=c
      real(rk) :: A(3,3)    = 0._rk
      real(rk) :: Ainv(3,3) = 0._rk      ! inversa para frac→cart
      ! --- control -------------------------------------------------
      logical  :: pbc(3) = [.true., .true., .true.]
      integer  :: dim    = 3            ! 2 (2-D) o 3 (3-D)
      logical  :: centered = .false.    ! frac en [0,1) o en [-0.5,0.5)
   contains
      procedure :: update     => cell_update
   end type Cell

   !======================== INTERFAZ PÚBLICA ========================
   public :: cell_from_lengths_angles
   public :: cell_from_vectors
   public :: cart_to_frac, frac_to_cart
   public :: wrap_frac_centered, wrap_cart, pbc_apply_point
   public :: min_image
   public :: pbc_ORTO        ! rutina legacy rápida ortorrómbica
   public :: cell_volume     ! determinante de la celda (volumen)


   ! alias legible
   interface pbc_apply_point
      module procedure wrap_cart
   end interface pbc_apply_point
   
contains
!====================================================================
subroutine cell_from_lengths_angles(c, a, b, c_, alpha, beta, gamma, dim, centered)
   type(Cell),      intent(out)           :: c
   real(rk),        intent(in)            :: a, b, c_
   real(rk),        intent(in)            :: alpha, beta, gamma
   integer,  optional, intent(in)         :: dim
   logical,  optional, intent(in)         :: centered
   !-----------------------------------------------------------------
   real(rk) :: alp, bet, gam, cosg, sing, cosb, sinb, cosa, sina
   real(rk) :: ax(3), bx(3), cx(3)        ! vectores

   ! longitudes y ángulos
   c%a_len     = a
   c%b_len     = b
   c%c_len     = c_
   c%alpha_deg = alpha
   c%beta_deg  = beta
   c%gamma_deg = gamma

   ! convierte a radianes
   alp = alpha * acos(-1._rk) / 180._rk
   bet = beta  * acos(-1._rk) / 180._rk
   gam = gamma * acos(-1._rk) / 180._rk
   cosa = cos(alp); sina = sin(alp)
   cosb = cos(bet); sinb = sin(bet)
   cosg = cos(gam); sing = sin(gam)

   ! vectores (convención estándar)
   ax = [ a, 0._rk, 0._rk ]
   bx = [ b*cosg, b*sing, 0._rk ]
   cx = [ c_*cosb,                                        &
          c_*(cosa - cosb*cosg)/sing,                     &
          c_*sqrt(1._rk - cosb**2 - ((cosa - cosb*cosg)/sing)**2) ]

   c%A(:,1) = ax
   c%A(:,2) = bx
   c%A(:,3) = cx

   c%dim      = merge(dim, 3, present(dim))
   c%centered = merge(centered, .false., present(centered))

   call c%update()        ! genera Ainv
end subroutine cell_from_lengths_angles
!====================================================================
subroutine cell_from_vectors(c, a_vec, b_vec, c_vec, dim, centered)
   type(Cell),      intent(out)           :: c
   real(rk),        intent(in)            :: a_vec(3), b_vec(3), c_vec(3)
   integer, optional, intent(in)          :: dim
   logical, optional, intent(in)          :: centered
   !-----------------------------------------------------------------
   c%A(:,1) = a_vec
   c%A(:,2) = b_vec
   c%A(:,3) = c_vec

   c%a_len = norm2(a_vec)
   c%b_len = norm2(b_vec)
   c%c_len = norm2(c_vec)

   c%alpha_deg = angle_deg(b_vec, c_vec)
   c%beta_deg  = angle_deg(a_vec, c_vec)
   c%gamma_deg = angle_deg(a_vec, b_vec)

   c%dim      = merge(dim, 3, present(dim))
   c%centered = merge(centered, .false., present(centered))

   call c%update()
contains
   pure real(rk) function norm2(v)
      real(rk), intent(in) :: v(:)
      norm2 = sqrt(sum(v*v))
   end function
   pure real(rk) function angle_deg(v1, v2)
      real(rk), intent(in) :: v1(:), v2(:)
      angle_deg = acos( dot_product(v1,v2) / (norm2(v1)*norm2(v2)) ) * 180._rk / acos(-1._rk)
   end function
end subroutine cell_from_vectors
!====================================================================
subroutine cell_update(c)
   class(Cell), intent(inout) :: c
   c%Ainv = inv3(c%A)
end subroutine cell_update
!====================================================================
pure function cart_to_frac(c, r_cart) result(s)
   type(Cell),  intent(in) :: c
   real(rk),    intent(in) :: r_cart(3)
   real(rk)                :: s(3)
   s = matmul(c%Ainv, r_cart)
   if (c%centered) s = s - nint(s)          ! lleva a [-0.5,0.5)
end function cart_to_frac
!====================================================================
pure function frac_to_cart(c, s) result(r_cart)
   type(Cell),  intent(in) :: c
   real(rk),    intent(in) :: s(3)
   real(rk)                :: r_cart(3)
   r_cart = matmul(c%A, s)
end function frac_to_cart
!====================================================================
pure function wrap_frac_centered(c, s) result(sw)
   type(Cell),  intent(in) :: c
   real(rk),    intent(in) :: s(3)
   real(rk)                :: sw(3)
   sw = s - nint(s)         ! [-0.5,0.5)
end function wrap_frac_centered
!====================================================================
pure function wrap_cart(c, r_cart) result(rw)
   type(Cell), intent(in)  :: c
   real(rk),   intent(in)  :: r_cart(3)
   real(rk)               :: rw(3)
   rw = frac_to_cart(c, wrap_frac_centered(c, cart_to_frac(c, r_cart)))
end function wrap_cart


!====================================================================

pure function min_image(c, r1, r2) result(dr_cart)
   type(Cell), intent(in) :: c
   real(rk),   intent(in) :: r1(3), r2(3)
   real(rk)              :: dr_cart(3)
   real(rk) :: s(3)

   ! 1) Vector diferencia en fraccionales (sin aplicar centered: no envolver aún)
   s = matmul(c%Ainv, r2 - r1)

   ! 2) Envoltura solo en ejes periódicos (slab: ejes no periódicos no se modifican)
   where (c%pbc)
      s = s - nint(s)                 ! => rango [-0.5,0.5)
   end where

   ! 3) Regreso a cartesianas
   dr_cart = frac_to_cart(c, s)
end function min_image

!====================================================================
!  Rutina rápida legacy ortorrómbica: lbox = [Lx,Ly,Lz]
!  ¡Solo válida si los vectores son ortogonales!
!====================================================================
pure function pbc_ORTO(dr, lbox) result(dr_min)
   real(rk), intent(in) :: dr(3), lbox(3)
   real(rk)             :: dr_min(3)
   dr_min = dr - lbox * nint(dr / lbox)
end function pbc_ORTO
!====================================================================
!  Volumen de la celda = det(A)
!====================================================================
pure real(rk) function cell_volume(c) result(vol)
   type(Cell), intent(in) :: c
   vol = c%A(1,1)*(c%A(2,2)*c%A(3,3) - c%A(2,3)*c%A(3,2)) &
       - c%A(1,2)*(c%A(2,1)*c%A(3,3) - c%A(2,3)*c%A(3,1)) &
       + c%A(1,3)*(c%A(2,1)*c%A(3,2) - c%A(2,2)*c%A(3,1))
end function cell_volume
!====================================================================
!  ────────── utilidades internas ──────────
!  Inversa 3×3 por cofactores                             (det ≠ 0)
!====================================================================
pure function inv3(M) result(invM)
   real(rk), intent(in) :: M(3,3)
   real(rk)             :: invM(3,3)
   real(rk) :: det
   invM(1,1) =  M(2,2)*M(3,3) - M(2,3)*M(3,2)
   invM(1,2) = -M(1,2)*M(3,3) + M(1,3)*M(3,2)
   invM(1,3) =  M(1,2)*M(2,3) - M(1,3)*M(2,2)

   invM(2,1) = -M(2,1)*M(3,3) + M(2,3)*M(3,1)
   invM(2,2) =  M(1,1)*M(3,3) - M(1,3)*M(3,1)
   invM(2,3) = -M(1,1)*M(2,3) + M(1,3)*M(2,1)

   invM(3,1) =  M(2,1)*M(3,2) - M(2,2)*M(3,1)
   invM(3,2) = -M(1,1)*M(3,2) + M(1,2)*M(3,1)
   invM(3,3) =  M(1,1)*M(2,2) - M(1,2)*M(2,1)

   det = M(1,1)*invM(1,1) + M(1,2)*invM(2,1) + M(1,3)*invM(3,1)
   invM = invM / det
end function inv3
!====================================================================
end module PBC_Mod
