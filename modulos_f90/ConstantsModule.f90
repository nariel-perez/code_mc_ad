MODULE PhysicalConstants
  USE PBC_Mod, only: rk
  IMPLICIT NONE

  private
  public  :: ComputeConstants, FCLEC, FACTORELEC
  REAL :: AK, AE0, FACTORELEC, FCLEC1AUX, FCLEC2AUX, FCLEC3AUX, FCLECAUX, FCLEC

CONTAINS
  SUBROUTINE ComputeConstants(SIGMETANO, SIGMA, EPS, AK_input, diel)
    IMPLICIT NONE
    REAL(rk), INTENT(IN) :: SIGMETANO, EPS, diel
    REAL, INTENT(IN) :: SIGMA, AK_input
    
    ! Compute physical constants
    AK = AK_input / 6.023E23
    AE0 = 8.85E-12   ! Permittivity in free space
    FACTORELEC = 96500.0 / 6.023E23
    FCLEC1AUX = SQRT(4.0 * 3.14159265 * AE0)
    FCLEC2AUX = SQRT(SIGMETANO / SIGMA * 1E-10)
    FCLEC3AUX = SQRT((EPS * AK) * diel)
    FCLECAUX = FCLEC1AUX * FCLEC2AUX * FCLEC3AUX
    FCLEC = 1.0 / FCLECAUX
    
    WRITE(*,*) ''
    WRITE(*,*) "-------------------------------------------------"
    WRITE(*,*) "------------ELECTRIC CONSTANTS-------------------"
    WRITE(*,*) ' '
    WRITE(*,*) 'e: ', FACTORELEC
    WRITE(*,*) ' FCLEC (to reduced units)', FCLEC
        
  END SUBROUTINE ComputeConstants
  
END MODULE PhysicalConstants
