!C------------------------------------------------------------------------
!C------------------------------------------------------------------------
!C------------------------------------------------------------------------
!C   SUBRUTINA ADPOTOUT
!C   CALCULA LA ENERGIA DE INTERACCION ENTRE LA SUPERFICIE Y EL FLUIDO
!C------------------------------------------------------------------------
!      
!C    *******************************************************************
!C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
!C    **                                                               **
!C    ** PRINCIPAL VARIABLES:                                          **
!C    **                                                               **
!C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
!C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
!C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
!C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
!C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
!C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
!C    ** REAL    SIGMA             LJ DIAMETER                         **
!C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
!C    **                                                               **
!C    ** USAGE:                                                        **
!C    **                                                               **
!C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
!C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
!C    ** RANGE CORRECTIONS IS INCLUDED.                                **
!C    *******************************************************************
!C    *******************************************************************
!----------------------------------------------------
!C     *******************************************************************

SUBROUTINE ADPOTOUT(IPULL, MOLKIND, DELTV)
    USE InputParams, ONLY: mat
    USE AdsorbateInput, ONLY: NATOMKIND, NATOM
    USE SimulationData, ONLY: RX, RY, RZ, UADS
    IMPLICIT NONE

    ! Argumentos de la subrutina
    INTEGER :: IPULL, MOLKIND
    REAL :: DELTV

    ! Variables locales
    INTEGER :: NMAX
    PARAMETER (NMAX = 5000)
    REAL :: RCUT, SIGMA
    REAL :: DELTV1, DELTW
    INTEGER :: NC
    REAL :: RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
    REAL :: RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
    REAL :: VLRC0, WLRC0
    INTEGER :: J, JIN, I1, I, K
    INTEGER :: IPOT

    ! Inicialización de acumuladores
    DELTV = 0.0
    DELTW = 0.0

    ! Bucle sobre todos los átomos de la molécula
    DO I1 = 1, NATOM(MOLKIND)
        RXI = RX(IPULL, I1, MOLKIND)
        RYI = RY(IPULL, I1, MOLKIND)
        RZI = RZ(IPULL, I1, MOLKIND)

        IPOT = NATOMKIND(I1, MOLKIND)

        ! Calcular índices de la celda
        I = INT(RXI * MAT)
        J = INT(RYI * MAT)
        K = INT(RZI * MAT)

        ! Calcular el potencial de adsorción
        DELTV1 = UADS(I, J, K, IPOT)
        DELTV = DELTV + DELTV1
    END DO

    ! Cambiar el signo de DELTV y DELTW para una eliminación
    DELTV = -DELTV
    DELTW = -DELTW

    RETURN
END SUBROUTINE ADPOTOUT
