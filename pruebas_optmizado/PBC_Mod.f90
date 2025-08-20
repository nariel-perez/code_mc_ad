!--------------------------------------------------------------------
! File: PBC_Mod.f90 (Refactorizado)
!
! Manejo genérico de celdas periódicas: construcción, conversiones
! cart ↔ frac, envolturas y vector mínima-imagen.
!--------------------------------------------------------------------
module PBC_Mod
    use, intrinsic :: iso_fortran_env, only : error_unit
    implicit none
    private ! Todo oculto; se exporta al final

    !============================= KIND ===============================
    integer, parameter, public :: rk = selected_real_kind(12, 99)

    !============================= TYPE ===============================
    type, public :: Cell
        ! --- datos geométricos ---------------------------------------
        ! Longitudes y ángulos definitorios de la celda
        real(rk) :: a_len = 0.0_rk, b_len = 0.0_rk, c_len = 0.0_rk
        real(rk) :: alpha_deg = 90.0_rk, beta_deg = 90.0_rk, gamma_deg = 90.0_rk
        
        ! Vectores de red en columnas: A(:,1)=a, A(:,2)=b, A(:,3)=c
        real(rk) :: A(3,3) = 0.0_rk
        
        ! Matriz inversa para la conversión de coordenadas
        real(rk) :: Ainv(3,3) = 0.0_rk
        
        ! --- control -------------------------------------------------
        ! Flag para condiciones de contorno periódicas en cada eje
        logical :: pbc(3) = [.true., .true., .true.]
        
        ! Dimensionalidad de la simulación
        integer :: dim = 3
        
        ! Flag para la envoltura de coordenadas: [0,1) o [-0.5,0.5)
        logical :: centered = .false.
    contains
        procedure :: update => cell_update
    end type Cell

    !======================== INTERFAZ PÚBLICA ========================
    public :: Cell, cell_from_lengths_angles, cell_from_vectors
    public :: cart_to_frac, frac_to_cart
    public :: wrap_frac, wrap_cart, pbc_apply_point
    public :: min_image
    public :: pbc_ORTO

    ! Alias para una interfaz más clara
    interface pbc_apply_point
        module procedure wrap_cart
    end interface pbc_apply_point
    
contains

!====================================================================
! Subrutinas para la construcción de la celda
!====================================================================

subroutine cell_from_lengths_angles(c, a, b, c_in, alpha, beta, gamma, dim, centered)
    ! Construye un objeto Cell a partir de las longitudes y ángulos de celda.
    !
    ! Args:
    ! c (out): El objeto Cell a construir.
    ! a, b, c_in (in): Longitudes de los vectores de la celda.
    ! alpha, beta, gamma (in): Ángulos de la celda en grados.
    ! dim (optional, in): Dimensionalidad de la simulación (2 o 3).
    ! centered (optional, in): Si las coordenadas fraccionales se centran.
    
    type(Cell), intent(out) :: c
    real(rk), intent(in) :: a, b, c_in
    real(rk), intent(in) :: alpha, beta, gamma
    integer, optional, intent(in) :: dim
    logical, optional, intent(in) :: centered
    
    real(rk) :: alp_rad, bet_rad, gam_rad
    real(rk) :: cos_gam, sin_gam, cos_bet, sin_bet, cos_alp, sin_alp
    
    ! Asignar parámetros de entrada
    c%a_len = a
    c%b_len = b
    c%c_len = c_in
    c%alpha_deg = alpha
    c%beta_deg = beta
    c%gamma_deg = gamma
    
    ! Convertir ángulos a radianes
    alp_rad = alpha * (acos(-1.0_rk) / 180.0_rk)
    bet_rad = beta * (acos(-1.0_rk) / 180.0_rk)
    gam_rad = gamma * (acos(-1.0_rk) / 180.0_rk)
    
    cos_alp = cos(alp_rad); sin_alp = sin(alp_rad)
    cos_bet = cos(bet_rad); sin_bet = sin(bet_rad)
    cos_gam = cos(gam_rad); sin_gam = sin(gam_rad)
    
    ! Construir la matriz de la celda (convención estándar)
    c%A(:,1) = [a, 0.0_rk, 0.0_rk]
    c%A(:,2) = [b*cos_gam, b*sin_gam, 0.0_rk]
    
    c%A(1,3) = c_in * cos_bet
    c%A(2,3) = c_in * (cos_alp - cos_bet*cos_gam) / sin_gam
    c%A(3,3) = c_in * sqrt(1.0_rk - cos_bet**2 - ((cos_alp - cos_bet*cos_gam)/sin_gam)**2)
    
    c%dim = merge(dim, 3, present(dim))
    c%centered = merge(centered, .false., present(centered))
    
    call c%update() ! Calcular la inversa
end subroutine cell_from_lengths_angles

subroutine cell_from_vectors(c, a_vec, b_vec, c_vec, dim, centered)
    ! Construye un objeto Cell a partir de los vectores de red.
    !
    ! Args:
    ! c (out): El objeto Cell a construir.
    ! a_vec, b_vec, c_vec (in): Vectores de la celda.
    ! dim (optional, in): Dimensionalidad.
    ! centered (optional, in): Si las coordenadas fraccionales se centran.
    
    type(Cell), intent(out) :: c
    real(rk), intent(in) :: a_vec(3), b_vec(3), c_vec(3)
    integer, optional, intent(in) :: dim
    logical, optional, intent(in) :: centered
    
    ! Funciones internas para calcular norma y ángulo
    pure real(rk) function norm2(v)
        real(rk), intent(in) :: v(:)
        norm2 = sqrt(sum(v*v))
    end function
    
    pure real(rk) function angle_deg(v1, v2)
        real(rk), intent(in) :: v1(:), v2(:)
        real(rk) :: dot, n1, n2
        dot = dot_product(v1, v2)
        n1 = norm2(v1)
        n2 = norm2(v2)
        ! Evitar divisiones por cero
        if (n1*n2 == 0.0_rk) then
            angle_deg = 0.0_rk
        else
            angle_deg = acos(dot / (n1*n2)) * 180.0_rk / acos(-1.0_rk)
        end if
    end function
    
    c%A(:,1) = a_vec
    c%A(:,2) = b_vec
    c%A(:,3) = c_vec
    
    c%a_len = norm2(a_vec)
    c%b_len = norm2(b_vec)
    c%c_len = norm2(c_vec)
    
    c%alpha_deg = angle_deg(b_vec, c_vec)
    c%beta_deg = angle_deg(a_vec, c_vec)
    c%gamma_deg = angle_deg(a_vec, b_vec)
    
    c%dim = merge(dim, 3, present(dim))
    c%centered = merge(centered, .false., present(centered))
    
    call c%update()
end subroutine cell_from_vectors

!====================================================================
! Métodos para la manipulación de la celda
!====================================================================

subroutine cell_update(c)
    ! Actualiza la matriz inversa Ainv del objeto Cell.
    ! Esta subrutina se llama automáticamente después de construir o modificar una celda.
    class(Cell), intent(inout) :: c
    c%Ainv = inv3(c%A)
end subroutine cell_update

pure function cart_to_frac(c, r_cart) result(s)
    ! Convierte un vector de coordenadas cartesianas a fraccionales.
    !
    ! Args:
    ! c (in): Objeto Cell.
    ! r_cart (in): Vector de coordenadas cartesianas (x, y, z).
    ! Returns:
    ! s: Vector de coordenadas fraccionales (a, b, c).
    
    type(Cell), intent(in) :: c
    real(rk), intent(in) :: r_cart(3)
    real(rk) :: s(3)
    
    s = matmul(c%Ainv, r_cart)
    if (c%centered) s = wrap_frac_centered(c, s)
end function cart_to_frac

pure function frac_to_cart(c, s) result(r_cart)
    ! Convierte un vector de coordenadas fraccionales a cartesianas.
    !
    ! Args:
    ! c (in): Objeto Cell.
    ! s (in): Vector de coordenadas fraccionales (a, b, c).
    ! Returns:
    ! r_cart: Vector de coordenadas cartesianas (x, y, z).
    
    type(Cell), intent(in) :: c
    real(rk), intent(in) :: s(3)
    real(rk) :: r_cart(3)
    
    r_cart = matmul(c%A, s)
end function frac_to_cart

pure function wrap_frac_centered(c, s) result(sw)
    ! Envuelve un vector de coordenadas fraccionales en el rango [-0.5, 0.5).
    ! Esta es la lógica central de envoltura.
    
    type(Cell), intent(in) :: c
    real(rk), intent(in) :: s(3)
    real(rk) :: sw(3)
    
    sw = s - nint(s)
end function wrap_frac_centered

pure function wrap_cart(c, r_cart) result(rw)
    ! Envuelve un vector de coordenadas cartesianas dentro de la celda.
    !
    ! Args:
    ! c (in): Objeto Cell.
    ! r_cart (in): Vector de coordenadas cartesianas.
    ! Returns:
    ! rw: Vector envuelto en la celda.
    
    type(Cell), intent(in) :: c
    real(rk), intent(in) :: r_cart(3)
    real(rk) :: rw(3)
    
    ! Conversión a fraccionales, envoltura y regreso a cartesianas
    rw = frac_to_cart(c, wrap_frac_centered(c, cart_to_frac(c, r_cart)))
end function wrap_cart

pure function min_image(c, r1, r2) result(dr_cart)
    ! Calcula el vector de distancia mínima entre dos puntos,
    ! considerando las condiciones de contorno periódicas.
    !
    ! Args:
    ! c (in): Objeto Cell con las propiedades de periodicidad.
    ! r1, r2 (in): Vectores de coordenadas cartesianas.
    ! Returns:
    ! dr_cart: Vector diferencia más corto (mínima imagen).
    
    type(Cell), intent(in) :: c
    real(rk), intent(in) :: r1(3), r2(3)
    real(rk) :: dr_cart(3), s(3)
    
    ! 1) Vector en coordenadas fraccionales
    s = cart_to_frac(c, r2 - r1)
    
    ! 2) Aplicar envoltura solo en ejes periódicos
    where (c%pbc)
        s = s - nint(s)
    end where
    
    ! 3) Regresar a coordenadas cartesianas
    dr_cart = frac_to_cart(c, s)
end function min_image

!====================================================================
! Rutinas de utilidad (pueden ser trasladadas a un módulo aparte)
!====================================================================

pure function pbc_ORTO(dr, lbox) result(dr_min)
    ! Rutina legacy rápida para celdas ortorrómbicas (solo si los vectores son ortogonales).
    !
    ! Args:
    ! dr (in): Vector diferencia entre dos puntos.
    ! lbox (in): Longitudes de la caja [Lx, Ly, Lz].
    ! Returns:
    ! dr_min: Vector de distancia mínima.
    
    real(rk), intent(in) :: dr(3), lbox(3)
    real(rk) :: dr_min(3)
    
    dr_min = dr - lbox * nint(dr / lbox)
end function pbc_ORTO

pure function inv3(M) result(invM)
    ! Calcula la inversa de una matriz 3x3 usando el método de cofactores.
    !
    ! Args:
    ! M (in): Matriz 3x3.
    ! Returns:
    ! invM: Matriz inversa.
    
    real(rk), intent(in) :: M(3,3)
    real(rk) :: invM(3,3)
    real(rk) :: det
    
    invM(1,1) = M(2,2)*M(3,3) - M(2,3)*M(3,2)
    invM(1,2) = -M(1,2)*M(3,3) + M(1,3)*M(3,2)
    invM(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
    invM(2,1) = -M(2,1)*M(3,3) + M(2,3)*M(3,1)
    invM(2,2) = M(1,1)*M(3,3) - M(1,3)*M(3,1)
    invM(2,3) = -M(1,1)*M(2,3) + M(1,3)*M(2,1)
    invM(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
    invM(3,2) = -M(1,1)*M(3,2) + M(1,2)*M(3,1)
    invM(3,3) = M(1,1)*M(2,2) - M(1,2)*M(2,1)
    
    det = M(1,1)*invM(1,1) + M(1,2)*invM(2,1) + M(1,3)*invM(3,1)
    if (abs(det) < epsilon(det)) then
        ! Podría ser una buena idea devolver un error en un contexto de programa
    end if
    
    invM = invM / det
end function inv3

end module PBC_Mod
