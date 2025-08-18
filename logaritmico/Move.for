C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        SUBROUTINE MOVE( TEMP, Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos)
c--------
        IMPLICIT NONE
	COMMON/BLOCK1/RX,RY,RZ,
     +NATOMKIND,EPSI,SIGM,Q,
     +RX0,RY0,RZ0,
     +RX1,RY1,RZ1,
     +RXC,RYC,RZC,
     +EPSAC,SGC,QAC,
     +UADS,acel,acelx,acely,acelz,
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
      REAL EXOLD,EYOLD,EZOLD
	REAL dELTVaDOUT2,DELTVOUT2

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
        REAL DELTAX,DELTAY,DELTAZ,DELTA
        REAL        RXBE(50),RYBE(50),RZBE(50)
	REAL        RXBE2(50),RYBE2(50),RZBE2(50)
        INTEGER I,IJ,J
        REAL VOLD,VGOLD,VAOLD
        REAL VANT,VGANT,VAANT
        REAL VNUEVA,VGNUEVA,VANUEVA
        REAL EX,EY,EZ,RR
        REAL        ANGUL
        PARAMETER (ANGUL=3.14159)
        REAL DX,DY,DZ
        REAL A11,A12,A13
        REAL A21,A22,A23
        REAL A31,A32,A33
        REAL DELTVADIN,DELTVIN
        REAL DELTVADOUT,DELTVOUT
        REAL RMIN
        REAL DELTCB
        REAL DELTANG
	real deltaxspc,deltayspc,deltazspc
	real rxnew,rynew,rznew
        
        INTEGER POSXIN,POSYIN,POSZIN
        LOGICAL TESTIGO,TESTIGO2
      real rxi,ryi,rzi
      
      real xmax,ymax,zmax
      !return
      
      xmax=acelx/acel
	!xmax=1
      ymax=acely/acel
	!ymax=1
      zmax=acelz/acel
	!zmax=1
C    *******************************************************************
	CREATE=.FALSE.
        GHOST  = .FALSE.
        BETA   = 1.0 / (TEMP)
      RMIN   = 0.75 * SIGMA

!-----------------------------------------------------------------
!     Elegir tipo de molécula Y RECUPERAR POSICIONES
!-----------------------------------------------------------------
	!return
        
        MOLKIND=RANF(DUMMY)*NMOLEC+1
        !molkind=1

        
        NTRIAL = N(MOLKIND) - 1       
        !write(*,*) 'MOVE N ',N(MOLKIND)
        
        IF ( NTRIAL .LT. 0 ) RETURN

C    ** PICK A RANDOM ELEMENT FROM THE ACTIVE PART OF LOCATE **
        b=RANF(DUMMY) 
        NLOC  = INT ( REAL ( NTRIAL ) * b ) + 1
          !nloc=1
        IPULL = LOCATE(NLOC,MOLKIND)
        !write(*,*) IPULL, ' IPULL1'
        !write(*,*)'POSICIONES INICIALES'
        DO I=1,NATOM(MOLKIND)
            RXBE(I)=RX(IPULL,I,MOLKIND)
            RYBE(I)=RY(IPULL,I,MOLKIND)
            RZBE(I)=RZ(IPULL,I,MOLKIND)
            
            RX1(I)=RXBE(I)
            RY1(I)=RYBE(I)
            RZ1(I)=RZBE(I)
            !write(*,*) rx1(i),ry1(i),rz1(i),i
            rxi=rx1(i)
            ryi=ry1(i)
            rzi=rz1(i)
        ENDDO
        !PAUSE
        !write(*,*) 'POSICIONES FINALES'
        EXOLD=0
        EYOLD=0
        EZOLD=0
        EXNEW=0
        EYNEW=0
        EZNEW=0
        IF(NATOM(MOLKIND).GT.1)THEN
        EXOLD=ANX(IPULL,MOLKIND)
        EYOLD=ANGY(IPULL,MOLKIND)
        EZOLD=ANZ(IPULL,MOLKIND)
        ENDIF

!-----------------------------------------------------------------------------
        !ENERGÍAS VIEJAS
        VOLD=V
        VGOLD=VG
        VAOLD=VA
!--------------------------------------------------------------
C    ** CALCULATE ENERGY CHANGE ON REMOVAL OF ATOM IPULL **
        CALL POTOUT ( IPULL, MOLKIND,  DELTVOUT)
        !write(*,*)DELTVOUT,'DELTVOUT'
      
	  CALL ADPOTOUT(IPULL, MOLKIND,DELTVADOUT )
      !write(*,*)DELTVADOUT,'DELTVADOUT'

          

!-------------------------------------------------------------
        !CALCULO ENERGÍA ANTES DEL PASO
        VANT=VOLD+DELTVOUT+DELTVADOUT
        VGANT=VGOLD+DELTVOUT
        VAANT=VAOLD+DELTVADOUT
        
!--------------------------------------------------------------
        !ELECCIÓN DEL TIPO DE MOVIMIENTO
        !IJ=1 ROTACIÓN
        !IJ=2 TRASLACIÓN
        
        IJ=INT(RANF(DUMMY)*2)+1
        IF(NATOM(MOLKIND).EQ.1) IJ=2
	!ij=1
        GOTO(10,20) IJ
		goto 20
!---------------------------------------------------------------------
!     MOLECULAR ROTATION
!---------------------------------------------------------------------
!     Elección del punto de una esfera
  10    IF(NATOM(MOLKIND).GT.1) THEN

       DO I=1,NATOM(MOLKIND)
            
            RXBE2(I)=RX0(I,MOLKIND)

            RYBE2(I)=RY0(I,MOLKIND)
            RZBE2(I)=RZ0(I,MOLKIND)
         !   !write(*,*) RXBE(I),RYBE(I),RZBE(I)
            
        ENDDO

	rxnew=rx1(2)
	rynew=ry1(2)
	rznew=rz1(2)

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
        
        DX=3.14159*EX
        DY=3.14159*EY
        DZ=3.14159*EZ
        
        A11=COS(DZ)*COS(DX)-SIN(DZ)*COS(DY)*SIN(DX)
        A12=SIN(DZ)*COS(DX)+COS(DZ)*COS(DY)*SIN(DX)
        A13=SIN(DY)*SIN(DX)
        A21=-COS(DZ)*SIN(DX)-SIN(DZ)*COS(DY)*COS(DX)
        A22=-SIN(DZ)*SIN(DX)+COS(DZ)*COS(DY)*COS(DX)
        A23=SIN(DY)*COS(DX)
        A31=SIN(DZ)*SIN(DY)
        A32=-COS(DZ)*SIN(DY)
        A33=COS(DY)
      DO I=1,NATOM(MOLKIND)
          RX1(I)=A11*RXBE2(I)+A12*RYBE2(I)+A13*RZBE2(I)
          RY1(I)=A21*RXBE2(I)+A22*RYBE2(I)+A23*RZBE2(I)
          RZ1(I)=A31*RXBE2(I)+A32*RYBE2(I)+A33*RZBE2(I)
      ENDDO

     	DO I=1,NATOM(MOLKIND)
          RX1(I)=RX1(I)+RXNEW
          RY1(I)=RY1(I)+RYNEW
          RZ1(I)=RZ1(I)+RZNEW
          
            RX1(I)=RX1(I)-BCX*xmax*ANINT(RX1(I)/xmax)
            RY1(I)=RY1(I)-BCY*ymax*ANINT(RY1(I)/ymax)
            RZ1(I)=RZ1(I)-BCZ*zmax*ANINT(RZ1(I)/zmax)
            
		!write(*,*) RX1(I),RY1(I),RZ1(I),i,' ** '
            IF (ABS(RX1(I)).GT.0.5) RETURN
       IF (ABS(RY1(I)).GT.0.5) RETURN
       IF (ABS(RZ1(I)).GT.0.5) RETURN
            
      END DO    
      ENDIF

      GOTO 30
!--------------------------------------------------------------
      !TRANSLACION
!--------------------------------------------------------------
        !CALCULO DESPLAZAMIENTO
   20     DELTA=0.01
        EXNEW=EXold
        EYNEW=EYold
        EZNEW=EZold

        DELTAX=(RANF(DUMMY)-0.5)*SIGMA*delta
        DELTAY=(RANF(DUMMY)-0.5)*SIGMA*delta
        DELTAZ=(RANF(DUMMY)-0.5)*SIGMA*delta
!      !write(*,*) 'DESPLAZAMIENTO'
      DO I=1,NATOM(MOLKIND)
          RX1(I)=RX1(I) +DELTAX
          RY1(I)=RY1(I) +DELTAY
          RZ1(I)=RZ1(I) +DELTAZ



           RX1(I)=RX1(I)-BCX*xmax*ANINT(RX1(I)/xmax)
           RY1(I)=RY1(I)-BCY*ymax*ANINT(RY1(I)/ymax)
           RZ1(I)=RZ1(I)-BCZ*zmax*ANINT(RZ1(I)/zmax)
        !write(*,*) rx1(i),ry1(i),rz1(i),i
          
       IF (ABS(RX1(I)).GT.0.5) RETURN
       IF (ABS(RY1(I)).GT.0.5) RETURN
       IF (ABS(RZ1(I)).GT.0.5) RETURN
      

      ENDDO
!--------------------------------------------------------------
   30          DO I=1,NATOM(MOLKIND)
            RX(IPULL,I,MOLKIND)=rx1(i)

            RY(IPULL,I,MOLKIND)=ry1(i)
            RZ(IPULL,I,MOLKIND)=rz1(i)
        	enddo
            ANX(IPULL,MOLKIND)=EXNEW
            ANGY(IPULL,MOLKIND)=EYNEW
            ANZ(IPULL,MOLKIND)=EZNEW
       
	CALL POTOUT ( IPULL, MOLKIND,  DELTVOUT2)
        
	  CALL ADPOTOUT(IPULL, MOLKIND,DELTVADOUT2 )
!--------------------------------------------------------------

        
!---------------------------------------------------------------
      !ENERGÍA NUEVA
      VNUEVA=VANT-DELTVout2-DELTVadout2
      VGNUEVA=VGANT-DELTVout2
      VANUEVA=VAANT-DELTVADout2
!----------------------------------------------------------------
      !CAMBIO EN LA ENERGÍA
      
      DELTCB=BETA*(VNUEVA-VOLD)
      !write(*,*) DELTCB,'****DELTCB'
      !PAUSE
      
!----------------------------------------------------------------
      IF (DELTCB.LT.75) THEN
          IF (DELTCB.le.0.0) THEN


              V=VNUEVA
              VG=VGNUEVA
              VA=VANUEVA
          ELSE IF ( EXP ( - DELTCB ) .GT. RANF(DUMMY) ) THEN
        
              V=VNUEVA
              VG=VGNUEVA
              VA=VANUEVA

          ELSE
                  DO I=1,NATOM(MOLKIND)
            RX(IPULL,I,MOLKIND)=rxbe(i)

            RY(IPULL,I,MOLKIND)=rybe(i)
            RZ(IPULL,I,MOLKIND)=rzbe(i)
	enddo
            ANX(IPULL,MOLKIND)=EXold
            ANGY(IPULL,MOLKIND)=EYold
            ANZ(IPULL,MOLKIND)=EZold
                  !write(*,*) 'NO MUEVO'
          V=VOLD
          VG=VGOLD
          VA=VAOLD
          ENDIF
      ELSE
                  DO I=1,NATOM(MOLKIND)
            RX(IPULL,I,MOLKIND)=rxbe(i)

            RY(IPULL,I,MOLKIND)=rybe(i)
            RZ(IPULL,I,MOLKIND)=rzbe(i)
	enddo
            ANX(IPULL,MOLKIND)=EXold
            ANGY(IPULL,MOLKIND)=EYold
            ANZ(IPULL,MOLKIND)=EZold

          V=VOLD
          VG=VGOLD
          VA=VAOLD

      ENDIF
	!write(*,*)'NMOLKIND=',N(MOLKIND)
      !!pause
      RETURN
        END
C------------------------------------------------------------------------------
