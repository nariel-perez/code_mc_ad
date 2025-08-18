C--------------------------------------------------------------------
      SUBROUTINE IN (TEMP,Z,SIGMA,EPS, RCUT, V,va,vg,W,CREATE,cr,
     +jpasos,canonicalmolecules)
      IMPLICIT NONE
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,ACEL, acelx,acely,acelz,
     +USS,FLAG,
     +BCX,BCY,BCZ,mat,
     +NMOLEC,N, NATOM,
     +ANX,ANGY,ANZ,EXNEW,EYNEW,EZNEW
     +	  /BLOCK2/LOCATE
     +	  /BLOCK3/nmaxi,nmin

C    *******************************************************************
C    ** ROUTINE TO ATTEMPT A TRIAL CREATION                           **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** REAL    TEMP            TEMPERATURE                           **
C    ** REAL    Z               ABSOLUTE ACTIVITY                     **
C    ** REAL    SIGMA           LENNARD-JONES DIAMETER                **
C    ** REAL    RCUT            CUT-OFF DISTANCE                      **
C    ** REAL    V               POTENTIAL ENERGY                      **
C    ** REAL    W               VIRIAL                                **
C    ** INTEGER N               NUMBER OF ATOMS BEFORE TRIAL CREATION **
C    ** LOGICAL CREATE          TRUE FOR A SUCCESSFUL CREATION        **
C    *******************************************************************
            REAL ANX(5000,10),ANGY(5000,10),ANZ(5000,10)
      REAL EXNEW,EYNEW,EZNEW
	integer canonicalmolecules
      REAL RXC(9000),RYC(9000),RZC(9000)
      REAL EPSAC(9000),SGC(9000),QAC(9000)
       INTEGER NATOMKIND (50,10)
        INTEGER     NMAX,NMIN(1000),NMAXI(1000)
	INTEGER NCONFMIN,NCONFMAX
        PARAMETER ( NMAX = 15000 )
        REAL        ANGUL
        PARAMETER (ANGUL=3.14159)
        REAL        TEMP, Z(10), SIGMA, RCUT
        REAL        RXBE(50),RYBE(50),RZBE(50)
        INTEGER     NC
        LOGICAL     CREATE
        INTEGER N(10)
        REAL        BETA, RXNEW, RYNEW, RZNEW
        REAL        RANF, DUMMY, RMIN
        INTEGER     NTRIAL
        LOGICAL     OVRLAP
        REAL XMAX,YMAX,ZMAX
        INTEGER MOLKIND
        INTEGER NMOLEC
        INTEGER I
        INTEGER NATOM(10)
        REAL RX0(50,10),RY0(50,10),RZ0(50,10),DX,DY,DZ
        REAL EX,EY,EZ,RR
        REAL COS,SIN
        REAL A11,A12,A13
        REAL A21,A22,A23
        REAL A31,A32,A33
        REAL RX1(50),RY1(50),RZ1(50)
        REAL DELTVA,V,VG,VA, DELTV, DELTW, DELTCB,W
        REAL CR
        REAL EPS
        REAL BCX,BCY,BCZ,acel,acelx,acely,acelz
        REAL RX(5000,50,10),RY(5000,50,10),RZ(5000,50,10)
        REAL EPSI(50),SIGM(50),Q(50)
        REAL UADS(-100:100,-100:100,-100:100,50)
      REAL USS(5000,50,50)
      LOGICAL FLAG(1000)
      INTEGER LOCATE(5000,10)
      INTEGER MAT,NTOTALB
	REAL ZTOTAL
      real P1,P2,P3
          INTEGER POSXIN,POSYIN,POSZIN,jpasos
!     ! WRITE(*,*)MAT, 'MAT IN'
!     !     write(*,*)ACEL, acelx,acely,acelz,' aceles******'
          EXNEW=0
          EYNEW=0
          EZNEW=0
        
		XMAX=acelx/acel
		YMAX=acely/acel
		ZMAX=acelz/acel


C    *******************************************************************

        CREATE = .FALSE.
        BETA   = 1.0 / (TEMP)
        RMIN   = 0.75 * SIGMA
        
!-----------------------------------------------------------------
!     Elegir tipo de molécula
!-----------------------------------------------------------------
        MOLKIND=INT(RANF(DUMMY)*nmolec)+1
!-----------------------------------------------------------------
        NTRIAL = N(MOLKIND) + 1
!-----------------------------------------------------------------
        IF ( NTRIAL .GE. NMAX ) return

C    ** GENERATE THE POSITION OF THE TRIAL ATOM CENTER OF MASS**



        RXNEW  = (RANF(DUMMY)-0.5)*XMAX
        RYNEW  = (RANF(DUMMY)-0.5)*YMAX
        RZNEW  = (RANF(DUMMY)-0.5)*ZMAX
        
	!write(*,*) 'generate'
	!write(*,*) rxnew, rynew,rznew,' POS'
	!write(*,*)'ntrial=',ntrial
	!pause
!------------------------------------------------------------------
!     GENERATE THE NEW POSITION
!------------------------------------------------------------------        
      
        DO I=1,NATOM(MOLKIND)
            
            RXBE(I)=RX0(I,MOLKIND)

            RYBE(I)=RY0(I,MOLKIND)
            RZBE(I)=RZ0(I,MOLKIND)
         !   write(*,*) RXBE(I),RYBE(I),RZBE(I)
            RX1(I)=RXBE(I)
            RY1(I)=RYBE(I)
            RZ1(I)=RZBE(I)
            
            
        ENDDO

!---------------------------------------------------------------------
!     MOLECULAR ROTATION
!---------------------------------------------------------------------
!     Elección del punto de una esfera
      IF(NATOM(MOLKIND).GT.1) THEN
        EX=(2.*RANF(DUMMY)-1.)
        EY=(2.*RANF(DUMMY)-1.)
        EZ=(2.*RANF(DUMMY)-1.)
        
        RR=SQRT(EX*EX+EY*EY+EZ*EZ)
        
        EX=EX/RR
        EY=EY/RR
        EZ=EZ/RR
        
        EXNEW=EX
        EYNEW=EY
        EZNEW=EZ
        
        DX=ANGUL*EX
        DY=ANGUL*EY
        DZ=ANGUL*EZ
        
        A11=COS(DZ)*COS(DX)-SIN(DZ)*COS(DY)*SIN(DX)
        A12=SIN(DZ)*COS(DX)+COS(DZ)*COS(DY)*SIN(DX)
        A13=SIN(DY)*SIN(DX)
        A21=-COS(DZ)*SIN(DX)-SIN(DZ)*COS(DY)*COS(DX)
        A22=-SIN(DZ)*SIN(DX)+COS(DZ)*COS(DY)*COS(DX)
        A23=SIN(DY)*COS(DX)
        A31=SIN(DZ)*SIN(DY)
        A32=-COS(DZ)*SIN(DY)
        A33=COS(DY)
       ! write(*,*) 'rotate'
      DO I=1,NATOM(MOLKIND)
          RX1(I)=A11*RXBE(I)+A12*RYBE(I)+A13*RZBE(I)
          RY1(I)=A21*RXBE(I)+A22*RYBE(I)+A23*RZBE(I)
          RZ1(I)=A31*RXBE(I)+A32*RYBE(I)+A33*RZBE(I)
	!write(*,*) RX1(I),RY1(I),RZ1(I)
   
      ENDDO
      ENDIF
	!write(*,*) ' IN 3'

      DO I=1,NATOM(MOLKIND)
          RX1(I)=RX1(I)+RXNEW
          RY1(I)=RY1(I)+RYNEW
          RZ1(I)=RZ1(I)+RZNEW
          
            RX1(I)=RX1(I)-BCX*xmax*ANINT(RX1(I)/xmax)
            RY1(I)=RY1(I)-BCY*ymax*ANINT(RY1(I)/ymax)
            RZ1(I)=RZ1(I)-BCZ*zmax*ANINT(RZ1(I)/zmax)
            
	!	!write(*,*) RX1(I),RY1(I),RZ1(I),' ** '
            IF (ABS(RX1(I)).GT.0.5) RETURN
       IF (ABS(RY1(I)).GT.0.5) RETURN
       IF (ABS(RZ1(I)).GT.0.5) RETURN
            
      ENDDO      
      !!write(*,*)molkind,' molkindin******3'

	CALL ADPOTIN (MOLKIND, DELTVA)
	!!write(*,*)ACEL, acelx,acely,acelz,' aceles******4'

	!write(*,*) deltva, ' deltva'
      if(deltva.gt.1000) return
      CALL POTIN ( MOLKIND, SIGMA,RCUT, RMIN, DELTV,
     :            OVRLAP )
	!write(*,*) deltv, ' deltv--5'

	!pause
      IF(OVRLAP) RETURN
C------------------------------------------------------------------------
C	CAMBIO EN EL LA INTERACCION SOLIDO-FLUIDO

C    ** CHECK FOR ACCEPTANCE **
      IF ( .NOT. OVRLAP ) THEN
          
          DELTCB = BETA*(DELTV+DELTVA)-LOG(z(molkind)/REAL(ntrial ))
	!!write(*,*) deltcb,' deltcb'
	!pause

          IF ( DELTCB .LT. 75.0 ) THEN
              IF ( DELTCB .LE. 0.0 ) THEN
                  CREATE = .TRUE.
                  CALL ADD ( MOLKIND )
                  V    = V + DELTV+DELTVA
                  VG=VG+deltV			
                  VA=VA+DELTVA
                  W    = W + DELTW
      !
                  N(MOLKIND)    = NTRIAL
		!pause
              ELSE IF ( EXP ( - DELTCB ) .GT. RANF ( DUMMY ) ) THEN
                  CREATE = .TRUE.
                  CALL ADD (  MOLKIND )
                  cr=cr+1
                  V    = V + DELTV + DELTVA
                  VG=VG+DELTV
                  VA=Va+DEltVA
                  W    = W + DELTW
                  N(MOLKIND)    = NTRIAL
			!!pause
              ENDIF
          ENDIF
      ENDIF
	!!write(*,*) ' IN 5'

      RETURN
      END
C------------------------------------------------------------------------------------
C------------------------------------------------------------------------------------
