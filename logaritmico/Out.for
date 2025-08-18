C------------------------------------------------------------------------------------
C------------------------------------------------------------------------------------

        SUBROUTINE OUT ( TEMP, Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos,
     +NMIN,NMAXI,canonicalmolecules)
c--------
        IMPLICIT NONE
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
C    ** ROUTINE TO ATTEMPT A TRIAL DESTRUCTION                        **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP         TEMPERATURE                              **
C    ** REAL    Z            ABSOLUTE ACTIVITY                        **
C    ** REAL    SIGMA        LENNARD-JONES DIAMETER                   **
C    ** REAL    RCUT         CUT-OFF DISTANCE                         **
C    ** REAL    V            POTENTIAL ENERGY                         **
C    ** REAL    W            VIRIAL                                   **
C    ** INTEGER N            NUMBER OF ATOMS BEFORE TRIAL DESTRUCTION **
C    ** LOGICAL GHOST        TRUE FOR A SUCCESSFUL DESTRUCTION        **
C    *******************************************************************
!----------------------------------------------------
            REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW

        INTEGER N(10),canonicalmolecules
       
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
      INTEGER NATOM(10),NMIN(1000),NMAXI(1000)
	INTEGER NCONFMIN,NCONFMAX
        INTEGER NTOTALB,I
      integer jpasos
C    *******************************************************************

        INTEGER     NMAX
        PARAMETER ( NMAX = 5000 )


        REAL        TEMP, Z(10), SIGMA, RCUT,EPS
        REAL      V, W, VA,VG

        INTEGER     NC
        LOGICAL     GHOST

        REAL        BETA, DELTV, DELTW, DELTDB, DELTVA
        REAL  RANF, DUMMY
        INTEGER     NTRIAL, NLOC, IPULL
        LOGICAL     OVRLAP,CREATE
        real sorteo, b
        real Ntotalmolec
	REAL ZTOTAL
        
        !RETURN
C    *******************************************************************
	CREATE=.FALSE.
        GHOST  = .FALSE.
        BETA   = 1.0 / (TEMP)
 !-----------------------------------------------------------------
!     Elegir tipo de molécula
!-----------------------------------------------------------------
        
        MOLKIND=RANF(DUMMY)*(NMOLEC)+1
        !molkind=1

        
        NTRIAL = N(MOLKIND) - 1       
        
        

!-----------------------------------------------------------------   
C    ** PICK A RANDOM ELEMENT FROM THE ACTIVE PART OF LOCATE **
        b=RANF(DUMMY) 
        NLOC  = INT ( REAL ( NTRIAL ) * b ) + 1
        IPULL = LOCATE(NLOC,MOLKIND)

C    ** CALCULATE ENERGY CHANGE ON REMOVAL OF ATOM IPULL **
!      WRITE(*,*) 'OUT', N(MOLKIND)
      
        CALL POTOUT ( IPULL, MOLKIND,  DELTV)
!        WRITE(*,*)DELTV,'DELTV'
      
	  CALL ADPOTOUT(IPULL, MOLKIND,DELTVA )
!      WRITE(*,*) DELTVA,'DELTVA'
C    ** CHECK FOR ACCEPTANCE **
       !Ntotalmolec=REAL(N(1)+N(2))
        DELTDB =BETA*(DELTV+DELTVA)-LOG(N(molkind)/Z(molkind))
!        WRITE(*,*) DELTDB,' DELTDB'
!------------------------------------------------------------------------
        IF ( DELTDB .LT. 75.0 ) THEN

           IF ( DELTDB .LT. 0.0 ) THEN

              GHOST = .TRUE.
              CALL REMOVE ( NLOC, IPULL,MOLKIND )

              V = V + DELTV+DELTVA
              VG=VG+DELTV
              VA=VA+DELTVA
              N(MOLKIND) = NTRIAL
              !write(*,*)'destruyo 1',n(molkind),deltdb,ipull
              !pause
           ELSE 
               sorteo=RANF ( DUMMY )
           IF ( EXP( -DELTDB ) .GT.sorteo  ) THEN
               


              GHOST = .TRUE.
              
              CALL REMOVE ( NLOC, IPULL,MOLKIND)
              

              V = V + DELTV+DELTVA
              
              VG=Vg+deltV
              VA=VA+DELTVA
              N(MOLKIND) = NTRIAL
              !write(*,*)'destruyo 2',n(molkind),sorteo,ipull
              !pause
           endif
           
           ENDIF

        ENDIF
	
	!!pause
      RETURN
        END
C------------------------------------------------------------------------------
