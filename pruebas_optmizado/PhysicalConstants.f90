!---------------------------------------------------------------------
! File: PhysicalConstants.f90 (Refactorizado)
!
! Módulo para la gestión y cálculo de constantes físicas.
! La inicialización se realiza de forma segura y controlada.
!---------------------------------------------------------------------
module PhysicalConstants
    use PBC_Mod, only: rk
    use, intrinsic :: iso_fortran_env, only : output_unit

    implicit none
    private
    
    !======================= VARIABLES PÚBLICAS =====================
    ! Se mantienen los nombres de variables originales
    public :: ComputeConstants, FCLEC, FACTORELEC
    
    ! Se declaran las constantes del módulo (ahora privadas)
    real(rk) :: AK, AE0, FACTORELEC, FCLEC

    !======================== CONSTANTES NOMBRADAS PRIVADAS ===========
    real(rk), parameter :: AVOGADRO = 6.02214076e23_rk
    real(rk), parameter :: PI = acos(-1.0_rk)
    
contains
    
    !------------------------------------------------------------------
    ! Procedimiento para el cálculo e inicialización de constantes
    !------------------------------------------------------------------
    subroutine ComputeConstants(SIGMETANO, SIGMA, EPS, AK_input, diel)
        ! Calcula y asigna los valores a las constantes del módulo.
        !
        ! Args:
        ! SIGMETANO, SIGMA, EPS (in): Parámetros de la simulación.
        ! AK_input (in): Constante de Boltzmann en unidades cal/mol K.
        ! diel (in): Constante dieléctrica del medio.
        
        implicit none
        real(rk), intent(in) :: SIGMETANO, SIGMA, EPS, AK_input, diel
        
        real(rk) :: FCLEC1AUX, FCLEC2AUX, FCLEC3AUX, FCLECAUX
        
        ! Calcular las constantes físicas
        AK = AK_input / AVOGADRO
        AE0 = 8.85e-12_rk
        FACTORELEC = 96500.0_rk / AVOGADRO
        
        ! Cálculos intermedios para FCLEC
        FCLEC1AUX = sqrt(4.0_rk * PI * AE0)
        FCLEC2AUX = sqrt(SIGMETANO / SIGMA * 1e-10_rk)
        FCLEC3AUX = sqrt((EPS * AK) * diel)
        FCLECAUX = FCLEC1AUX * FCLEC2AUX * FCLEC3AUX
        FCLEC = 1.0_rk / FCLECAUX
        
        ! Impresión de resultados
        write(output_unit, '(A)') ''
        write(output_unit, '(A)') "-------------------------------------------------"
        write(output_unit, '(A)') "------------ELECTRIC CONSTANTS-------------------"
        write(output_unit, '(A)') ' '
        write(output_unit, '(A,ES15.8)') 'e: ', FACTORELEC
        write(output_unit, '(A,ES15.8)') 'FCLEC (to reduced units): ', FCLEC
        
    end subroutine ComputeConstants
    
end module PhysicalConstants
