!---------------------------------------------------------------------
! File: GeomUtils.f90 (Refactorizado)
!
! Módulo de utilidades geométricas para simulaciones con condiciones
! de contorno periódicas (PBC). Incluye funciones para calcular
! la matriz métrica y la distancia de mínima imagen.
!---------------------------------------------------------------------
module GeomUtils
    use PBC_Mod, only : rk, Cell_t => Cell
    
    implicit none
    private
    public :: cell_to_metric, r2_min_image_frac, wrap_by_pbc

contains
    
    !====================================================================
    ! Funciones de conversión y cálculo de métricas
    !====================================================================
    pure subroutine cell_to_metric(c, G, Ainv, pbcx, pbcy, pbcz)
        ! Obtiene la matriz métrica (G), la inversa de la matriz de la celda
        ! (Ainv) y los indicadores de periodicidad (pbcx, pbcy, pbcz) de un objeto Cell.
        !
        ! Args:
        ! c (in): El objeto Cell con las propiedades de la celda.
        ! G (out): La matriz métrica G = A^T * A.
        ! Ainv (out): La matriz inversa de la celda.
        ! pbcx, pbcy, pbcz (out): Indicadores booleanos de periodicidad para cada eje.
        
        type(Cell_t), intent(in) :: c
        real(rk), dimension(3,3), intent(out) :: G, Ainv
        logical, intent(out) :: pbcx, pbcy, pbcz
        
        Ainv = c%Ainv
        G = matmul(transpose(c%A), c%A)
        pbcx = c%pbc(1)
        pbcy = c%pbc(2)
        pbcz = c%pbc(3)
        
    end subroutine cell_to_metric
    
    pure real(rk) function r2_min_image_frac(G, s1, s2, s3) result(r2)
        ! Calcula el cuadrado de la distancia de mínima imagen entre dos
        ! puntos en coordenadas fraccionales.
        !
        ! Args:
        ! G (in): La matriz métrica de la celda.
        ! s1, s2, s3 (in): Coordenadas fraccionales de un vector diferencia.
        ! Returns:
        ! r2: El cuadrado de la distancia de mínima imagen.
        
        real(rk), dimension(3,3), intent(in) :: G
        real(rk), intent(in) :: s1, s2, s3
        real(rk), dimension(3) :: s_vec
        
        s_vec = (/ s1, s2, s3 /)
        r2 = matmul(transpose(s_vec), matmul(G, s_vec))
        
    end function r2_min_image_frac
    
    pure subroutine wrap_by_pbc(s1, s2, s3, pbcx, pbcy, pbcz)
        ! Envuelve las coordenadas fraccionales de un vector en el rango
        ! [-0.5, 0.5) para los ejes periódicos.
        !
        ! Args:
        ! s1, s2, s3 (inout): Coordenadas fraccionales del vector a envolver.
        ! pbcx, pbcy, pbcz (in): Indicadores booleanos de periodicidad.
        
        real(rk), intent(inout) :: s1, s2, s3
        logical, intent(in) :: pbcx, pbcy, pbcz
        
        if (pbcx) s1 = s1 - nint(s1)
        if (pbcy) s2 = s2 - nint(s2)
        if (pbcz) s3 = s3 - nint(s3)
        
    end subroutine wrap_by_pbc
    
end module GeomUtils
