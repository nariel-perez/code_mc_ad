!---------------------------------------------------------------------
! File: RotationModule.f90 (Refactorizado)
!
! Módulo para la gestión de rotaciones 3D. Utiliza tablas precalculadas
! para una evaluación rápida de funciones trigonométricas.
!---------------------------------------------------------------------
module RotationModule
    use PBC_mod, only: rk
    
    implicit none
    private
    
    !======================== CONSTANTES PÚBLICAS ====================
    ! `PI` y el tamaño de la tabla son públicos para la inicialización
    ! desde el programa principal si fuera necesario.
    integer, parameter, public :: TABLE_SIZE = 4096 ! Incrementado para mayor precisión
    real(rk), parameter, public :: PI = acos(-1.0_rk)
    
    !======================== VARIABLES DEL MÓDULO ===================
    ! La tabla de valores precalculados ahora es privada, gestionada
    ! por las subrutinas del módulo.
    real(rk), dimension(TABLE_SIZE) :: TABLE_SIN, TABLE_COS
    logical :: is_initialized = .false.
    
    !======================== INTERFAZ PÚBLICA =======================
    public :: init_rotation_tables, get_rotation_matrix

contains

!====================================================================
! Subrutinas de inicialización y utilidades
!====================================================================
subroutine init_rotation_tables()
    ! Inicializa las tablas precalculadas de seno y coseno.
    ! Este procedimiento es idempotente: solo se ejecutará una vez.
    
    implicit none
    integer :: i
    real(rk) :: angle, step_size
    
    if (is_initialized) return
    
    step_size = (2.0_rk * PI) / (TABLE_SIZE - 1)
    
    do i = 1, TABLE_SIZE
        ! El rango de ángulos es [-PI, PI]
        angle = -PI + (real(i - 1) * step_size)
        TABLE_SIN(i) = sin(angle)
        TABLE_COS(i) = cos(angle)
    end do
    
    is_initialized = .true.
    
end subroutine init_rotation_tables

!--------------------------------------------------------------------
pure function get_trig_value(angle, table) result(value)
    ! Busca un valor trigonométrico en una tabla con interpolación lineal
    ! para mejorar la precisión.
    !
    ! Args:
    ! angle (in): Ángulo en radianes en el rango [-PI, PI].
    ! table (in): La tabla a consultar (TABLE_SIN o TABLE_COS).
    ! Returns:
    ! value: El valor trigonométrico interpolado.
    
    real(rk), intent(in) :: angle
    real(rk), dimension(:), intent(in) :: table
    real(rk) :: value
    integer :: index
    real(rk) :: x1, x2, y1, y2, frac
    
    ! Normalizar el ángulo al rango de la tabla [-PI, PI]
    angle = angle - 2.0_rk * PI * nint(angle / (2.0_rk * PI))

    ! Calcular el índice base para la interpolación
    index = int((angle + PI) / ((2.0_rk * PI) / (TABLE_SIZE - 1))) + 1

    ! Clamp the index to prevent out-of-bounds access
    index = max(1, min(TABLE_SIZE - 1, index))

    ! Puntos para la interpolación
    x1 = -PI + (real(index - 1) * (2.0_rk * PI) / (TABLE_SIZE - 1))
    x2 = -PI + (real(index) * (2.0_rk * PI) / (TABLE_SIZE - 1))
    y1 = table(index)
    y2 = table(index + 1)
    
    ! Interpolación lineal
    frac = (angle - x1) / (x2 - x1)
    value = y1 + frac * (y2 - y1)
    
end function get_trig_value

!====================================================================
! Subrutina principal para obtener la matriz de rotación
!====================================================================
subroutine get_rotation_matrix(angle_x, angle_y, angle_z, R)
    ! Calcula la matriz de rotación 3x3 a partir de los ángulos de Euler
    ! (Z-Y'-X'' o Z-Y-X). Utiliza tablas precalculadas con interpolación.
    !
    ! Args:
    ! angle_x, angle_y, angle_z (in): Ángulos de rotación en radianes.
    ! R (out): Matriz de rotación 3x3 resultante.
    
    implicit none
    real(rk), intent(in) :: angle_x, angle_y, angle_z
    real(rk), dimension(3,3), intent(out) :: R
    
    real(rk) :: sen_dx, cos_dx, sen_dy, cos_dy, sen_dz, cos_dz
    
    ! Asegurar que las tablas estén inicializadas
    if (.not. is_initialized) then
        call init_rotation_tables()
    end if
    
    ! Obtener valores trigonométricos de la tabla con interpolación
    sen_dx = get_trig_value(angle_x, TABLE_SIN)
    cos_dx = get_trig_value(angle_x, TABLE_COS)
    sen_dy = get_trig_value(angle_y, TABLE_SIN)
    cos_dy = get_trig_value(angle_y, TABLE_COS)
    sen_dz = get_trig_value(angle_z, TABLE_SIN)
    cos_dz = get_trig_value(angle_z, TABLE_COS)
    
    ! Calcular la matriz de rotación Z-Y'-X'' (extrinsic) o Z-Y-X (intrinsic)
    ! Esta es la convención estándar en química computacional.
    R(1,1) = cos_dz * cos_dy
    R(1,2) = cos_dz * sen_dy * sen_dx - sen_dz * cos_dx
    R(1,3) = cos_dz * sen_dy * cos_dx + sen_dz * sen_dx

    R(2,1) = sen_dz * cos_dy
    R(2,2) = sen_dz * sen_dy * sen_dx + cos_dz * cos_dx
    R(2,3) = sen_dz * sen_dy * cos_dx - cos_dz * sen_dx

    R(3,1) = -sen_dy
    R(3,2) = cos_dy * sen_dx
    R(3,3) = cos_dy * cos_dx
    
end subroutine get_rotation_matrix

end module RotationModule
