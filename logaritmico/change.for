C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        SUBROUTINE change(TEMP,Z,SIGMA,EPS, RCUT,V,va,vg,W,GHOST,jpasos)
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
C    ** REAL    TEMP         TEMPERATURE                              **33
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

        INTEGER N(10)
       
        INTEGER MOLKIND1, molkind2
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
        REAL RXOR,RYOR,RZOR
        REAL RXNEW,RYNEW,RZNEW

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
        REAL B1,B2,B3
        INTEGER POSXIN,POSYIN,POSZIN
        LOGICAL TESTIGO,TESTIGO2
      real rxi,ryi,rzi
      
      real xmax,ymax,zmax
      
      
      !return
	if(nmolec.lt.2) return
	
	
	xmax=acelx/acel
	ymax=acely/acel
	zmax=acelz/acel
	
	
C    *******************************************************************
	CREATE=.FALSE.
        GHOST  = .FALSE.
        BETA   = 1.0 / (TEMP)
      RMIN   = 0.75 * SIGMA

!-----------------------------------------------------------------
!     Elegir tipo de molécula Y RECUPERAR POSICIONES
!-----------------------------------------------------------------
!        write(*,*)nmolec, n(1:3), '******* 1 '
        MOLKIND1=RANF(DUMMY)*NMOLEC+1
        !molkind=1

        
        NTRIAL = N(MOLKIND1) - 1       
        !WRITE(*,*) 'MOVE N ',N(MOLKIND)
        
        IF ( NTRIAL .LT. 0 ) RETURN

C    ** PICK A RANDOM ELEMENT FROM THE ACTIVE PART OF LOCATE **
        b=RANF(DUMMY) 
        NLOC  = INT ( REAL ( NTRIAL ) * b ) + 1
          !nloc=1
        IPULL = LOCATE(NLOC,MOLKIND1)
!        WRITE(*,*) IPULL, ' IPULL1'
!        write(*,*)'POSICIONES INICIALES'
        DO I=1,NATOM(MOLKIND1)
            RXBE(I)=RX(IPULL,I,MOLKIND1)
            RYBE(I)=RY(IPULL,I,MOLKIND1)
            RZBE(I)=RZ(IPULL,I,MOLKIND1)
            
            RX1(I)=RXBE(I)
            RY1(I)=RYBE(I)
            RZ1(I)=RZBE(I)
!            write(*,*) rx1(i),ry1(i),rz1(i),i
            rxi=rx1(i)
            ryi=ry1(i)
            rzi=rz1(i)
        ENDDO
        !----------------------------------------------------
        ! Determinación del origen
        RXOR=RXBE(1)-RX0(1,MOLKIND1)
        RYOR=RYBE(1)-RY0(1,MOLKIND1)
        RZOR=RZBE(1)-RZ0(1,MOLKIND1)

        !PAUSE
!        WRITE(*,*) 'POSICIONES FINALES'
        EXOLD=0
        EYOLD=0
        EZOLD=0
        EXNEW=0
        EYNEW=0
        EZNEW=0
        IF(NATOM(MOLKIND1).GT.1)THEN
        EXOLD=ANX(IPULL,MOLKIND1)
        EYOLD=ANGY(IPULL,MOLKIND1)
        EZOLD=ANZ(IPULL,MOLKIND1)
        ENDIF

!-----------------------------------------------------------------------------
        !ENERGÍAS VIEJAS
        VOLD=V
        VGOLD=VG
        VAOLD=VA
!	write(*,*)nmolec, n(1:3), '******* 2 '
!--------------------------------------------------------------
C    ** CALCULATE ENERGY CHANGE ON REMOVAL OF ATOM IPULL **
        CALL POTOUT ( IPULL, MOLKIND1,  DELTVOUT)
!        WRITE(*,*)DELTVOUT,'DELTVOUT'
      
	  CALL ADPOTOUT(IPULL, MOLKIND1,DELTVADOUT )
!      WRITE(*,*)DELTVADOUT,'DELTVADOUT'

          

!-------------------------------------------------------------
        !CALCULO ENERGÍA ANTES DEL PASO
        VANT=VOLD+DELTVOUT+DELTVADOUT
        VGANT=VGOLD+DELTVOUT
        VAANT=VAOLD+DELTVADOUT
!--------------------------------------------------------------
        !Creaci´n de la nueva molécula
        
  97    MOLKIND2=INT(RANF(DUMMY)*NMOLEC)+1
!  97	MOLKIND2=INT(RANF(DUMMY)*2)+1
!	IF (MOLKIND2.EQ.2) MOLKIND2=-1
!        MOLKIND2=molkind1+MOLKIND2!
!	if (molkind2.gt.nmolec.or.molkind2.lt.1) then
!	MOLKIND2=INT(RANF(DUMMY)*NMOLEC)+1
!	endif

!	write(*,*)molkind2,nmolec, n(1:3), '******* 3 '

        if(molkind1.eq.molkind2) goto 97
        
!-----------------------------------------------------------------

C    ** GENERATE THE POSITION OF THE TRIAL ATOM CENTER OF MASS**



        RXNEW  = RXOR
        RYNEW  = RYOR
        RZNEW  = RZOR
        
!------------------------------------------------------------------
!     GENERATE THE NEW POSITION
!------------------------------------------------------------------        
      
        DO I=1,NATOM(MOLKIND2)
            
            RXBE2(I)=RX0(I,MOLKIND2)
            RYBE2(I)=RY0(I,MOLKIND2)
            RZBE2(I)=RZ0(I,MOLKIND2)
            
            RX1(I)=RXBE2(I)
            RY1(I)=RYBE2(I)
            RZ1(I)=RZBE2(I)
            
            
        ENDDO

!---------------------------------------------------------------------
!     MOLECULAR ROTATION
!---------------------------------------------------------------------
!     Elección del punto de una esfera
      IF(NATOM(MOLKIND2).GT.1) THEN
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
        
      DO I=1,NATOM(MOLKIND2)
          RX1(I)=A11*RXBE2(I)+A12*RYBE2(I)+A13*RZBE2(I)
          RY1(I)=A21*RXBE2(I)+A22*RYBE2(I)+A23*RZBE2(I)
          RZ1(I)=A31*RXBE2(I)+A32*RYBE2(I)+A33*RZBE2(I)
   
      ENDDO
      ENDIF
!	write(*,*)nmolec, n(1:3), '******* 4 '

      DO I=1,NATOM(MOLKIND2)
          RX1(I)=RX1(I)+RXNEW
          RY1(I)=RY1(I)+RYNEW
          RZ1(I)=RZ1(I)+RZNEW
!          WRITE(*,*)I,RX1(I)*39.5,RY1(I)*39.5,RZ1(I)*39.5,ntrial,'NTR'
          
            RX1(I)=RX1(I)-BCX*xmax*ANINT(RX1(I)/xmax)
            RY1(I)=RY1(I)-BCY*ymax*ANINT(RY1(I)/ymax)
            RZ1(I)=RZ1(I)-BCZ*zmax*ANINT(RZ1(I)/zmax)
            
            IF (ABS(RX1(I)).GT.0.5) RETURN
       IF (ABS(RY1(I)).GT.0.5) RETURN
       IF (ABS(RZ1(I)).GT.0.5) RETURN
            
      ENDDO      

!--------------------------------------------------------------
 
  30  CALL ADPOTIN (MOLKIND2, DELTVADIN)
!      WRITE(*,*)DELTVADIN, 'DELTVADIN'
      
      !PAUSE
      IF (DELTVADIN.GT.100) RETURN
!--------------------------------------------------------------
        !REMUEVO LA MOLECULA
      
        CALL REMOVE (NLOC,IPULL,MOLKIND1)
        N(MOLKIND1)=N(MOLKIND1)-1
!--------------------------------------------------------------
!	write(*,*)nmolec, n(1:3), '******* 5 '

        
        !CALCULO NUEVA ENERGIA
        CALL POTIN ( MOLKIND2, SIGMA,RCUT, RMIN, DELTVIN,
     :            OVRLAP )
        IF(OVRLAP) THEN
            !WRITE(*,*) 'OVRLAP'
            !PAUSE
        DO I=1,NATOM(MOLKIND1)
            RX1(I)=RXBE(I)
            RY1(I)=RYBE(I)
            RZ1(I)=RZBE(I)
            
        ENDDO
        EXNEW=EXOLD
        EYNEW=EYOLD
        EZNEW=EZOLD
        CALL ADD(MOLKIND1)
        N(MOLKIND1)=N(MOLKIND1)+1
        RETURN
        ENDIF
!      WRITE(*,*)DELTVIN,' DELTVIN'
!---------------------------------------------------------------
      !ENERGÍA NUEVA
      VNUEVA=VANT+DELTVADIN+DELTVIN
      VGNUEVA=VGANT+DELTVIN
      VANUEVA=VAANT+DELTVADIN
!----------------------------------------------------------------
      !CAMBIO EN LA ENERGÍA
!      WRITE(*,*) DELTCB,'DELTCB'
      !PAUSE
      B1=Z(MOLKIND2)*REAL(N(MOLKIND1)+1)
      B2=Z(MOLKIND1)*REAL(N(MOLKIND2)+1)
      B3=LOG(B1/B2)
      !DELTCB=BETA*(VNUEVA-VOLD)
!---------------------------------------------------------------------------
!****************************************************************************
	!B3=0

      DELTCB = BETA*(VNUEVA-VOLD)-B3
!	write(*,*)nmolec, n(1:3), '******* 6 '

!----------------------------------------------------------------
      IF (DELTCB.LT.75) THEN
          IF (DELTCB.LE.0.0) THEN

              CALL ADD(MOLKIND2)
              V=VNUEVA
              VG=VGNUEVA
              VA=VANUEVA
              N(MOLKIND2)=N(MOLKIND2)+1
!		write(*,*)' muevo 1'
          ELSE IF ( EXP ( - DELTCB ) .GT. RANF(DUMMY) ) THEN
!	write(*,*)nmolec, n(1:3), '******* 6 bis'

          CALL ADD(MOLKIND2)
 !         	write(*,*)nmolec, n(1:3), '******* 6 tris'

              V=VNUEVA
!	write(*,*)nmolec, n(1:3), '******* 6 1'
              VG=VGNUEVA
!	write(*,*)nmolec, n(1:3), '******* 6 2'

              VA=VANUEVA
!	write(*,*)nmolec, n(1:3), '******* 6 3'
!	write(*,*)N(MOLKIND2), molkind2,' molkind2'
	
              N(MOLKIND2)=N(MOLKIND2)+1
!	write(*,*)N(MOLKIND2), molkind2,' molkind2'
!	write(*,*)nmolec, n(1:3), '******* 6 4'

 !             WRITE(*,*) 'MUEVO2'
              !PAUSE
!	write(*,*)nmolec, n(1:3), '******* 6 4-tris'
          ELSE
 !                           WRITE(*,*) 'NO MUEVO'
              !PAUSE

              DO I=1,NATOM(MOLKIND1)
                  RX1(I)=RXBE(I)
                  RY1(I)=RYBE(I)
                  RZ1(I)=RZBE(I)
              ENDDO
        EXNEW=EXOLD
        EYNEW=EYOLD
        EZNEW=EZOLD

          CALL ADD(MOLKIND1)
          V=VOLD
          VG=VGOLD
          VA=VAOLD
          
          N(MOLKIND1)=N(MOLKIND1)+1
          ENDIF
      ELSE
  !                          WRITE(*,*) 'NO MUEVO2'
              !PAUSE
!	write(*,*)nmolec, n(1:3), '******* 7 '
		!pause

          DO I=1,NATOM(MOLKIND1)
                  RX1(I)=RXBE(I)
                  RY1(I)=RYBE(I)
                  RZ1(I)=RZBE(I)
          ENDDO
        EXNEW=EXOLD
        EYNEW=EYOLD
        EZNEW=EZOLD
!	write(*,*)nmolec, n(1:3), '******* 7bis '

          CALL ADD(MOLKIND1)                    
          V=VOLD
          VG=VGOLD
          VA=VAOLD

          N(MOLKIND1)=N(MOLKIND1)+1
      ENDIF
 !     pause
!	write(*,*)nmolec, n(1:3), '******* 8 '
 
       RETURN
        END
C------------------------------------------------------------------------------
