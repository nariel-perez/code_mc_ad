!C-------------------------------------------------------------------------
!C-------------------------------------------------------------------------
!C------------------------------------------------------------------------
!C   SUBRUTINA ADPOTIN
!C   CALCULA LA ENERGIA DE INTERACCION ENTRE LA SUPERFICIE Y EL FLUIDO
!    C------------------------------------------------------------------------
!      
!C    *******************************************************************
!C    ** RETURNS THE POTENTIAL ENERGY CHANGE WITH SURFACE ON ADDING AN ATOM.    !    **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE ADDITION **
!C    ** REAL    RXI,RYI,RZI       THE COORDINATES OF THE ADDED ATOM   **
!C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
!C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
!C    ** REAL    SIGMA             LJ DIAMETER                         **
!C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
!C    ** REAL    RMIN              MINIMUM ALLOWED APPROACH OF ATOMS   **
!C    ** LOGICAL OVRLAP            TRUE FOR SUBSTANTIAL ATOM OVERLAP   **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!C    ** DURING A TRIAL ADDITION OF AN ATOM TO THE FLUID. THE LONG     **
!C    ** RANGE CORRECTIONS IS INCLUDED.                                **
!C    *******************************************************************


SUBROUTINE ADPOTIN(MOLKIND, DELTV)
    USE InputParams, ONLY: mat
    USE AdsorbateInput, ONLY: NATOM, NATOMKIND
    USE SimulationData, ONLY: RX1, RY1, RZ1,  UADS
    
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER, INTENT(IN) :: MOLKIND
    REAL, INTENT(OUT)   :: DELTV

    ! Variables locales
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    REAL :: RCUT, RMIN, SIGMA, RXI, RYI, RZI
    REAL :: DELTV1, DELTW
    INTEGER :: NC, NTRIAL
    LOGICAL :: OVRLAP
    REAL :: RCUTSQ, RMINSQ, SIGSQ, SR2, SR6, SR3, SR9
    REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
    REAL :: VLRC0, WLRC0
    INTEGER :: J, JIN
    INTEGER :: I1
    INTEGER :: IPOT
    INTEGER :: I, K

    ! Inicialización de variables
    DELTV = 0.0
    DELTW = 0.0
    DELTV1 = 0.0

    ! Bucle sobre los átomos de la molécula
    DO I1 = 1, NATOM(MOLKIND)
        RXI = RX1(I1)
        RYI = RY1(I1)
        RZI = RZ1(I1)

        ! Obtener el tipo de átomo
        IPOT = NATOMKIND(I1, MOLKIND)

        ! Calcular índices de la celda
        I = INT(RXI * MAT)
        J = INT(RYI * MAT)
        K = INT(RZI * MAT)

        ! Calcular el potencial de adsorción
        DELTV1 = UADS(I, J, K, IPOT)

        ! Verificar si el potencial es demasiado grande
        IF (DELTV1 > 100) THEN
            DELTV = 1E10
            RETURN
        END IF

        ! Acumular el potencial total
        DELTV = DELTV + DELTV1
    END DO

    RETURN
END SUBROUTINE ADPOTIN
