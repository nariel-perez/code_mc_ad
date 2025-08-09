!******************************************************************************
!* FICHE F.13.  THE HEART OF A CONSTANT MU VT MONTE CARLO PROGRAM             **
!* This FORTRAN code is intended to illustrate points made in the text.       **
!* To our knowledge it works correctly.  However it is the responsibility of  **
!* the user to test it, if it is to be used in a research application.        **
!******************************************************************************

!C    *******************************************************************
!C    ** ATTEMPTED CREATIONS AND DESTRUCTIONS IN GRAND CANONICAL MC.   **
!C    **                                                               **
!C    ** THESE ROUTINES ALLOW FOR A TRIAL DESTRUCTION OR CREATION IN A **
!C    ** GRAND CANONICAL MONTE CARLO PROGRAM.                          **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                   NUMBER OF ATOMS BEFORE THE TRIAL  **
!C    ** INTEGER NTRIAL              NUMBER OF ATOMS DURING THE TRIAL  **
!C    ** INTEGER NMAX                MAXIMUM NUMBER OF ATOMS ALLOWED   **
!C    ** INTEGER LOCATE(NMAX)        ARRAY OF ACTIVE ATOM INDICES      **
!C    ** REAL    RXNEW,RYNEW,RZNEW   POSITION FOR ADDITION OF ATOM     **
!C    ** REAL    RX(NMAX) ETC.       POSITIONS OF CURRENT ATOMS        **
!C    ** REAL    V                   POTENTIAL ENERGY + LRC            **
!C    ** REAL    W                   VIRIAL + LRC                      **
!C    ** REAL    DELTV               CHANGE IN ENERGY                  **
!C    ** REAL    DELTW               CHANGE IN VIRIAL                  **
!C    ** REAL    TEMP                REDUCED TEMPERATURE               **
!C    ** REAL    Z                   ABSOLUTE ACTIVITY COEFFICIENT     **
!C    ** REAL    SIGMA               LENNARD JONES DIAMETER            **
!C    ** REAL    RCUT                REDUCED CUTOFF DISTANCE           **
!C    ** REAL    RMIN                REDUCED MINIMUM SEPARATION        **
!C    ** LOGICAL OVRLAP              TRUE FOR SUBSTANTIAL ATOM OVERLAP **
!C    ** LOGICAL CREATE              TRUE FOR AN ACCEPTED CREATION     **
!C    ** LOGICAL GHOST               TRUE FOR AN ACCEPTED DESTRUCTION  **
!C    **                                                               **
!C    ** ROUTINES SUPPLIED:                                            **
!C    **                                                               **
!C    ** SUBROUTINE IN ( TEMP, Z, SIGMA, RCUT, N, V, W, CREATE )       **
!C    **    PERFORMS A TRIAL CREATION                                  **
!C    ** SUBROUTINE OUT ( TEMP, Z, SIGMA, RCUT, N, V, W, GHOST )       **
!C    **    PERFORMS A TRIAL DESTRUCTION                               **
!C    ** SUBROUTINE POTIN ( RXNEW, RYNEW, RZNEW, N, SIGMA, RCUT, RMIN, **
!C    ** :                  DELTV, DELTW, OVRLAP )                     **
!C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON CREATION         **
!C    ** SUBROUTINE POTOUT ( IPULL, N, SIGMA, RCUT, DELTV, DELTW )     **
!C    **    CALCULATES THE POTENTIAL ENERGY CHANGE ON DESTRUCTION      **
!C    **                                                               **
!C    ** ROUTINES REFERENCED:                                          **
!C    **                                                               **
!C    **    REAL FUNCTION RANF ( DUMMY )  (GIVEN IN F.11)              **
!C    **    RETURNS A UNIFORM RANDOM VARIATE ON ZERO TO ONE            **
!C    ** SUBROUTINE ADD ( RXNEW, RYNEW, RZNEW, N )                     **
!C    **    UPDATES LOCATE AFTER ADDITION (GIVEN IN F.14)              **
!C    ** SUBROUTINE REMOVE ( NLOC, IPULL, N )                          **
!C    **    UPDATES LOCATE AFTER REMOVAL (GIVEN IN F.14)               **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** ROUTINES IN AND OUT SHOULD BE CALLED WITH EQUAL PROBABILITY   **
!C    ** IN A GRAND CANONICAL MONTE CARLO SIMULATION. IF A TRIAL       **
!C    ** CREATION IS ACCEPTED THEN CREATE IS SET TO TRUE. IF A TRIAL   **
!C    ** DESTRUCTION IS ACCEPTED THEN GHOST IS SET TO TRUE. THE        **
!C    ** ROUTINES ARE WRITTEN FOR LENNARD-JONES ATOMS. THE BOX IS OF   **
!C    ** UNIT LENGTH, ALL DISTANCES ARE SCALED TO THE BOX LENGTH.      **
!C    ** TRIAL INPUTS WHICH RESULT IN A SEPARATION OF LESS THAN        **
!C    ** 0.5*SIGMA ARE REJECTED. THE LONG-RANGE CORRECTIONS ARE        **
!C    ** INCLUDED IN V AND W. ALL ACCUMULATORS ARE UPDATED IN THE MAIN **
!C    ** PART OF THE PROGRAM WHICH IS NOT GIVEN HERE.                  **
!C    *******************************************************************

program main
  use InputParams
  use simulationdata
  use AdsorbateInput
  use physicalconstants
  use Estructuramodule
  use rotationmodule
  
   implicit none
   
   ! ================================================================
   ! ============== Declaración de variables ========================
   ! ================================================================

   real :: p_ratio, SIGMA, TEMP, PRED, RCUT
   integer  :: i, NC,auxmat
   real, allocatable :: p_vals(:), Z(:)
   real :: XMAX, YMAX, ZMAX, VOL
   
   real, parameter   :: AK_input = 8.31   ! constante de los gases en j/mol·K
   




   
   ! ================================================================
   ! =================== Inicio del programa =========================
   ! ================================================================


   ! ----------------------------------------------------------------
   ! LECTURA DE ENTRADA
   ! ----------------------------------------------------------------

   call read_input('input.txt')  

   call print_params()
   
   
   auxmat = int(mat/2)
   allocate(uads(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, 50))
   
   !---------------------------------------
   !  pasos logaritmicos
   !--------------------------------------
   
   if (.not. allocated(p_vals)) allocate(p_vals(isot +1 ))
   p_vals = 0 
   
   p_ratio = (dp/p)**(1.0d0/ real(isot -1, kind = 8))

   p_vals(1) = p
   do  i = 2, isot
      p_vals(i) = p_vals(i-1)*p_ratio
   end do
   
   
   ! ----------------------------------------------------------------
   ! TRANSFORMACIÓN A UNIDADES REDUCIDAS
   ! ----------------------------------------------------------------
   SIGMA = sigmetano / acel
   TEMP  = T / eps
   P     = P ! * 1333.22
   PRED  = P * SIGMA**3 / eps
   RCUT  = 10*SIGMA
   
   XMAX = acelx / acel
   YMAX = acely / acel
   ZMAX = acelz / acel
   VOL  = XMAX * YMAX * ZMAX
   
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '--------------REDUCED UNITS ---------------------'
   write(*,'(A, F10.4)') 'TEMPERATURE: ', TEMP
   write(*,'(A, ES12.5)') 'PRESSURE: ', PRED
   
   write(*,'(A, F10.4)') 'SIGMA:', SIGMA
   write(*,'(A, ES12.4)') 'VOLUME:', VOL

   
   ! ----------------------------------------------------------------
   ! CÁLCULO DE CONSTANTES ELÉCTRICAS Y TABLAS DE ROTACIÓN
   ! ----------------------------------------------------------------
   call computeconstants(sigmetano, SIGMA, eps, AK_input, diel)
   call initrotationtables()

   ! ----------------------------------------------------------------
   ! LECTURA DE ESTRUCTURAS MOLECULARES
   ! ----------------------------------------------------------------
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '----------------ADSORBATES-----------------------'

   call read_adsorbates('MOLEC.DAT',sigma,  1.0e-7)
   
   if (.not. allocated(Z))    allocate(Z(NMOLEC))
   if (.not. allocated(N))    allocate(N(NMOLEC))
   
   if (.not. allocated(EPSI)) allocate(EPSI(maxAtoms))
   if (.not. allocated(SIGM)) allocate(SIGM(maxAtoms))
   if (.not. allocated(Q))    allocate(Q(maxAtoms))
   
   if (.not. allocated(RX1))  allocate(RX1(maxAtoms))
   if (.not. allocated(RY1))  allocate(RY1(maxAtoms))
   if (.not. allocated(RZ1))  allocate(RZ1(maxAtoms))
   if (.not. allocated(RX))   allocate(RX(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RY))   allocate(RY(5000, maxAtoms, NMOLEC))
   if (.not. allocated(RZ))   allocate(RZ(5000, maxAtoms, NMOLEC))
   if (.not. allocated(USS))  allocate(USS(5000, maxAtoms, maxAtoms))
   
   auxmat = int(NCELLMAT/2)
   
   if (.not. allocated(CNF)) then
      allocate(CNF(-auxmat:auxmat, -auxmat:auxmat, -auxmat:auxmat, NMOLEC, maxAtoms))
   end if
   
   
   
   ! ----------------------------------------------------------------
   ! LLAMADA A SUBRUTINAS DE POTENCIALES
   ! ----------------------------------------------------------------
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '------------------SURFACE------------------------'
   
   open(unit=50, file='SALIDAACTIVADO-100.TXT')
   open(unit=97, file='PERFILES.TXT')
   call estructura(eps, nam, sigma, sigmetano, NC, diel)
   
   write(*,*)
   write(*,*) '-------------------------------------------------'
   write(*,*) '------------POTENTIALS--------------------'
   write(*,*)
   
   call potencialff(eps, sigma, sigmetano, NC, RCUT, diel)
   call potencial(eps, sigma, sigmetano, NC, RCUT, diel)
   
   
end program main
