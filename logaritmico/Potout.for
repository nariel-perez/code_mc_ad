C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------


        SUBROUTINE POTOUT ( IPULL, MOLKIND,DELTV)
        IMPLICIT NONE
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,    acel,acelx,acely,acelz,
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
      REAL BCX,BCY,BCZ ,acel,acelx,acely,acelz
      INTEGER MAT
      integer NMOLEC
      INTEGER NATOM(10)
        

C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )

        REAL        RCUT, SIGMA, DELTV, DELTW
        INTEGER     IPULL

        REAL        RCUTSQ, SIGSQ, SR2, SR6, SR3, SR9, RXI, RYI, RZI
        REAL        RXIJ, RYIJ, RZIJ, RIJSQ, VIJ, WIJ, SIGCUB, PI
        REAL        VLRC0, WLRC0
        INTEGER     J, JIN

        PARAMETER ( PI = 3.14159265 )
        INTEGER I1
        INTEGER IPOT,JPOT
        INTEGER I
        INTEGER K
        REAL RIJ
        integer IDIST
        real xmax,ymax,zmax
        
        
        
        
        xmax=acelx/acel
        ymax=acely/acel
        zmax=acelz/acel
C     ******************************************************************
C    ** ZERO ACCUMULATORS **

C    ** ZERO ACCUMULATORS **

        DELTV  = 0.0
        DELTW  = 0.0

C    ** LOOP OVER ALL ATOMS  EXCEPT THE SAME KIND**
      !write(*,*) ipull, 'ipull potput'
      DO I1=1,NATOM(MOLKIND)
          
          RXI=RX(IPULL,I1,MOLKIND)
          RYI=RY(IPULL,I1,MOLKIND)
          RZI=RZ(IPULL,I1,MOLKIND)
          
          !write(*,*) ipull, rxi,ryi,rzi
          !write(*,*) 'pot out'
          
          IPOT=NATOMKIND(I1,MOLKIND)
          !write(*,*) ipot,' ipot'
          !write(*,*) nmolec, 'nmolec'
          !pause
          DO  I = 1,NMOLEC
             ! write(*,*)molkind,i,' molkind 1'
              IF(I.NE.MOLKIND) THEN
                  DO J=1,N(I)
                  !     ** PICK ACTIVE ATOMS FROM THE ARRAY LOCATE **

                  JIN   = LOCATE(J,I)
                  
                  DO K=1,NATOM(I)
                      JPOT=NATOMKIND(K,I)
                      RXIJ  = RXI - RX(JIN,K,I)
                      RYIJ  = RYI - RY(JIN,K,I)
                      RZIJ  = RZI - RZ(JIN,K,I)
                      
                      RXIJ  = RXIJ - BCX*xmax*ANINT ( RXIJ/xmax )
                      RYIJ  = RYIJ - BCY*ymax*ANINT ( RYIJ/ymax )
                      RZIJ  = RZIJ - BCZ*zmax*ANINT ( RZIJ/zmax )
                      
                      RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                      
                      RIJ=SQRT(RIJSQ)
                      
                      IDIST=RIJ*1000+1
                      
                      VIJ   = USS(IDIST,JPOT,IPOT)
                      DELTV = DELTV + VIJ
                      DELTW = DELTW + WIJ
                  ENDDO
                  ENDDO
                  else
                      
                  DO J=1,N(I)
!     ** PICK ACTIVE ATOMS FROM THE ARRAY LOCATE **
                  
                  JIN   = LOCATE(J,I)
                  !write(*,*)jin,ipull, ' jin ipull'
                  !write(*,*)molkind,i,' molkind'
                  IF(JIN.NE.IPULL) THEN
                  DO K=1,NATOM(I)
                      JPOT=NATOMKIND(K,I)
                      RXIJ  = RXI - RX(JIN,K,I)
                      RYIJ  = RYI - RY(JIN,K,I)
                      RZIJ  = RZI - RZ(JIN,K,I)
                      
                      RXIJ  = RXIJ - BCX*xmax*ANINT ( RXIJ/xmax )
                      RYIJ  = RYIJ - BCY*ymax*ANINT ( RYIJ/ymax )
                      RZIJ  = RZIJ - BCZ*zmax*ANINT ( RZIJ/zmax )
                      
                      RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
                      
                      RIJ=SQRT(RIJSQ)
                      
                      IDIST=RIJ*1000+1
                      
                      VIJ   = USS(IDIST,JPOT,IPOT)
                      !write(*,*)vij,idist,' vij out'
                      !write(*,*) ipot,jpot, 'ipot,jpot'
                      !pause
                      DELTV = DELTV + VIJ
                      DELTW = DELTW + WIJ
                  ENDDO
                  ENDIF
                  enddo
                  ENDIF
                  
                  
              ENDDO
          ENDDO
        DELTV = DELTV

C    ** CHANGE SIGN OF DELTV AND DELTW FOR A REMOVAL **

        DELTV = - DELTV
        DELTW = - DELTW
        RETURN
        END
C------------------------------------------------------------------------------
