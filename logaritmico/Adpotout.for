C------------------------------------------------------------------------
C------------------------------------------------------------------------
C------------------------------------------------------------------------
C   SUBRUTINA ADPOTOUT
C   CALCULA LA ENERGIA DE INTERACCION ENTRE LA SUPERFICIE Y EL FLUIDO
C------------------------------------------------------------------------
      SUBROUTINE ADPOTOUT (IPULL, MOLKIND, DELTV)
      implicit none
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS, acel,acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
C    *******************************************************************
C    ** RETURNS THE POTENTIAL ENERGY CHANGE WHEN AN ATOM IS REMOVED.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N                 THE NUMBER OF ATOMS BEFORE REMOVAL  **
C    ** INTEGER IPULL             THE ATOM TO BE REMOVED              **
C    ** INTEGER LOCATE(NMAX)      ARRAY OF ACTIVE ATOM INDICES        **
C    ** REAL    RX(NMAX) ETC.     THE ATOM POSITIONS                  **
C    ** REAL    DELTV             THE CHANCE IN POTENTIAL             **
C    ** REAL    DELTW             THE CHANGE IN VIRIAL                **
C    ** REAL    SIGMA             LJ DIAMETER                         **
C    ** REAL    RCUT              CUTOFF DISTANCE FOR POTENTIAL       **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** THIS SUBROUTINE IS USED TO CALCULATE THE CHANGE OF ENERGY     **
C    ** DURING A TRIAL DELETION OF AN ATOM FROM THE FLUID. THE LONG   **
C    ** RANGE CORRECTIONS IS INCLUDED.                                **
C    *******************************************************************
C    *******************************************************************
!----------------------------------------------------
            REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW

        INTEGER N(10)
       
        INTEGER MOLKIND
        REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
        
        INTEGER NATOMKIND(50,10)
        
      REAL UADS(-100:100,-100:100,-100:100,50)
      REAL USS(5000,50,50)
      LOGICAL FLAG(1000)
      INTEGER LOCATE(5000,10)
      REAL RXC(9000),RYC(9000),RZC(9000)
      REAL EPSAC(9000),SGC(9000),QAC(9000)
      REAL RX1(50),RY1(50),RZ1(50)
      REAL EPSI(50),SIGM(50),Q(50)
      REAL RX0(50,10),RY0(50,10),RZ0(50,10)
      REAL BCX,BCY,BCZ,acel,acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
        

C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA 
        REAL DELTV1, DELTV, DELTW
        INTEGER     NC, IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN, I1, I,K
        INTEGER IPOT


C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0
        DO I1=1,NATOM(MOLKIND)
            RXI=RX(IPULL,I1,MOLKIND)
            RYI=RY(IPULL,I1,MOLKIND)
            RZI=RZ(IPULL,I1,MOLKIND)
            
            IPOT=NATOMKIND(I1,MOLKIND)
            !write(*,*) rxi,ryi,rzi,' RX en ADPOT', ipot,' ipot'
            I=INT(RXI*mat)
            J=INT(RYI*mat)
            K=INT(RZI*mat)
            !write(*,*)i,j,k, 'ijk'
            DELTV1=UADS(I,J,K,IPOT)
            DELTV=DELTV+DELTV1
        ENDDO

        DELTV = - DELTV
        DELTW = - DELTW
        RETURN
        END

C-----------------------------------------------------------------------
